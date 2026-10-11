# Input data index

`input_index.csv` lists 71 inputs from the study's analysis collection.
Their roles vary across analyses. The index distinguishes 29 files with
official file or release locators, 20 derived inputs, eight entries linked
to study DOIs, and 14 entries for which an exact public location has not been
established. These categories are mutually exclusive. Some derived inputs also
have a DOI linking to their source study, so the number of rows with a URL
differs from the number of directly downloadable files.

For derived inputs, the index gives preparation notes and generating code
where available. A study DOI identifies the upstream publication, not a
downloadable copy of an intermediate file. The listed analysis paths describe
the expected input layout; files not distributed here must be obtained or
prepared separately.

Five ENCODE CTCF and DNase files were checked against the official file API,
including assembly, file size and MD5. GENCODE v19, the hg38-to-hg19 chain and
100-mer mappability references were checked against official release pages.
Their listed SHA-256 values identify the analysis inputs; these three files
were not downloaded again for packaging.

Replication timing uses a consensus of Repli-seq tracks from 15 cell lines
and the HeLa-S3 and HepG2 component tracks. The index gives their official
filenames and hashes.

Source data retain the access conditions and licenses of the originating
studies. Checksums identify files; they do not confer redistribution rights.
