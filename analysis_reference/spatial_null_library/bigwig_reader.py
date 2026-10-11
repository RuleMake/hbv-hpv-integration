# Read-only BigWig definition, original script main removed.
from __future__ import annotations
import struct
import zlib
from pathlib import Path
import numpy as np
BIGWIG_MAGIC = 0x888FFC26
BPT_MAGIC = 0x78CA8C91
RTREE_MAGIC = 0x2468ACE0
class BigWig:
    def __init__(self, path: str | Path):
        self.f = open(path, "rb")
        magic = struct.unpack("<I", self.f.read(4))[0]
        if magic != BIGWIG_MAGIC:
            raise ValueError(f"not a little-endian bigWig (magic {magic:#x})")
        (self.version, self.zoom_levels, self.chrom_tree_off, self.full_data_off,
         self.full_index_off, self.field_count, self.defined_field_count,
         self.autosql_off, self.total_summary_off, self.uncompress_buf_size,
         _reserved) = struct.unpack("<HHQQQHHQQIQ", self.f.read(60))

        self.zooms = []
        for _ in range(self.zoom_levels):
            red, _res, doff, ioff = struct.unpack("<IIQQ", self.f.read(24))
            self.zooms.append(dict(reduction=red, data_off=doff, index_off=ioff))

        self.f.seek(self.total_summary_off)
        vc, mn, mx, sd, ss = struct.unpack("<Qdddd", self.f.read(40))
        self.total_summary = dict(valid_count=vc, min=mn, max=mx, sum_data=sd, sum_sq=ss)
        self.chroms = self._read_chrom_tree()

    # ---------------------------------------------------------------- chrom B+ tree
    def _read_chrom_tree(self):
        self.f.seek(self.chrom_tree_off)
        magic, _bs, key_size, _vs, _n, _r = struct.unpack("<IIIIQQ", self.f.read(32))
        if magic != BPT_MAGIC:
            raise ValueError("bad chromosome B+ tree magic")
        out = {}

        def node():
            is_leaf, _pad, count = struct.unpack("<BBH", self.f.read(4))
            if is_leaf:
                for _ in range(count):
                    key = self.f.read(key_size).rstrip(b"\x00").decode()
                    cid, csize = struct.unpack("<II", self.f.read(8))
                    out[key] = dict(id=cid, size=csize)
            else:
                offs = []
                for _ in range(count):
                    self.f.read(key_size)
                    offs.append(struct.unpack("<Q", self.f.read(8))[0])
                for o in offs:
                    self.f.seek(o)
                    node()

        node()
        return out

    # ---------------------------------------------------------------- R-tree leaves
    def _leaves(self, index_off: int):
        self.f.seek(index_off)
        magic = struct.unpack("<I", self.f.read(4))[0]
        if magic != RTREE_MAGIC:
            raise ValueError("bad R-tree magic")
        self.f.read(44)  # blockSize, itemCount, start/end chrom+base, endFileOffset, itemsPerSlot, reserved
        root = self.f.tell()
        found = []

        def walk(off):
            self.f.seek(off)
            is_leaf, _pad, count = struct.unpack("<BBH", self.f.read(4))
            if is_leaf:
                raw = self.f.read(32 * count)
                for i in range(count):
                    _sc, _sb, _ec, _eb, doff, dsize = struct.unpack_from("<IIIIQQ", raw, 32 * i)
                    found.append((doff, dsize))
            else:
                raw = self.f.read(24 * count)
                kids = [struct.unpack_from("<IIIIQ", raw, 24 * i)[4] for i in range(count)]
                for k in kids:
                    walk(k)

        walk(root)
        return found

    def _zoom_records(self, zoom):
        """Yield (chrom_id, start, end, valid_count, sum_data) for one zoom level."""
        for doff, dsize in self._leaves(zoom["index_off"]):
            self.f.seek(doff)
            buf = self.f.read(dsize)
            if self.uncompress_buf_size > 0:
                buf = zlib.decompress(buf)
            n = len(buf) // 32
            for i in range(n):
                cid, st, en, vc, _mn, _mx, sd, _ss = struct.unpack_from("<IIIIffff", buf, 32 * i)
                yield cid, st, en, vc, sd

    # ---------------------------------------------------------------- public
    def pick_zoom(self, target: int):
        ok = [z for z in self.zooms if z["reduction"] <= target]
        if not ok:
            raise ValueError(f"no zoom level at or below {target} bp")
        return max(ok, key=lambda z: z["reduction"])

    def binned_means(self, bin_size: int = 50_000, target_reduction: int | None = None):
        """Mean value per fixed bin, aggregated from a zoom level. Returns a dict
        {chrom_name: (bin_starts, means)} plus the global mean actually observed."""
        zoom = self.pick_zoom(target_reduction or bin_size)
        id2name = {v["id"]: k for k, v in self.chroms.items()}
        acc = {k: dict(sd=np.zeros(v["size"] // bin_size + 2),
                       vc=np.zeros(v["size"] // bin_size + 2))
               for k, v in self.chroms.items()}
        tot_sd = tot_vc = 0.0
        for cid, st, en, vc, sd in self._zoom_records(zoom):
            name = id2name.get(cid)
            if name is None or vc == 0:
                continue
            a = acc[name]
            b = st // bin_size
            if b < len(a["sd"]):
                a["sd"][b] += sd
                a["vc"][b] += vc
            tot_sd += sd
            tot_vc += vc
        out = {}
        for name, a in acc.items():
            with np.errstate(invalid="ignore", divide="ignore"):
                m = np.where(a["vc"] > 0, a["sd"] / np.maximum(a["vc"], 1), np.nan)
            out[name] = (np.arange(len(m)) * bin_size, m)
        return out, (tot_sd / tot_vc if tot_vc else np.nan), zoom["reduction"]

    def validate(self, tol: float = 0.02):
        """Zoom-derived genome-wide mean must reproduce the header totalSummary mean."""
        hdr = self.total_summary["sum_data"] / self.total_summary["valid_count"]
        _, zoom_mean, red = self.binned_means(50_000)
        ok = abs(zoom_mean - hdr) / abs(hdr) <= tol
        return dict(header_mean=hdr, zoom_mean=zoom_mean, reduction=red,
                    rel_diff=abs(zoom_mean - hdr) / abs(hdr), passed=bool(ok))

    def close(self):
        self.f.close()