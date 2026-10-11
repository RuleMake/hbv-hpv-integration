"""Check donor-proportion statistics against constructed examples and saved results."""
import unittest
import numpy as np
import pandas as pd
import fig5c_donor_duplication_comparison as m


class EventUnitTests(unittest.TestCase):
    def test_two_endpoints_boundary_and_absent_chromosome(self):
        ev = pd.DataFrame({"chrom": ["chr1"]*3+["chr2"], "left_pos0": [90,100100,100101,100], "right_pos0": [110,200000,200001,120]})
        co = pd.DataFrame({"chrom": ["chr1"], "pos0": [100]})
        got = m.score_events(ev,co)
        self.assertEqual(got.near.tolist(), [True,True,False,False])
        self.assertEqual(int(got.near.sum()), 2)
        self.assertEqual(int((got.left_hit & got.right_hit).sum()), 1)

    def test_duplicate_coordinates_and_endpoint_swap(self):
        ev = pd.DataFrame({"chrom": ["chr1","chr1"], "left_pos0": [0,250000], "right_pos0": [100000,300000]})
        co = pd.DataFrame({"chrom": ["chr1"], "pos0": [1]})
        original = m.score_events(ev,co).near.tolist()
        self.assertEqual(original, m.score_events(ev,pd.concat([co,co],ignore_index=True)).near.tolist())
        swapped = ev.rename(columns={"left_pos0":"right_pos0","right_pos0":"left_pos0"})
        self.assertEqual(original,m.score_events(swapped,co).near.tolist())


class StatisticalTests(unittest.TestCase):
    def test_constant_equal_groups(self):
        groups = [(np.array([.5,.5]), np.array([.5,.5,.5]))]
        ci_d, ci_r, bad = m.bootstrap_effects(groups,np.array([1.]),np.random.default_rng(1),500)
        np.testing.assert_array_equal(ci_d,[0,0])
        np.testing.assert_array_equal(ci_r,[1,1])
        self.assertEqual(bad,0)
        self.assertEqual(m.permutation_test(groups,np.array([1.]),np.random.default_rng(1),100), (1.,1.))

    def test_exact_small_permutation_reference(self):
        # Six labeled individuals, choose three: 20 allocations; 2 maximally separated.
        groups = [(np.zeros(3),np.ones(3))]
        p_t,p_diff = m.permutation_test(groups,np.array([1.]),np.random.default_rng(2),19999)
        self.assertAlmostEqual(p_t,.1,delta=.015)
        self.assertAlmostEqual(p_diff,.1,delta=.015)

    def test_standardization_is_not_pooled_ratio(self):
        groups = [(np.array([0.,1.]),np.array([0.,0.])),(np.array([.5,.5]),np.array([1.,1.]))]
        a,b,d,r,t = m.weighted_point(groups,np.array([.25,.75]))
        self.assertEqual(a,.5)
        self.assertEqual(b,.75)
        self.assertEqual(d,-.25)
        self.assertAlmostEqual(r,2/3)

    def test_zero_bootstrap_denominator_is_exposed(self):
        groups=[(np.ones(2),np.zeros(2))]
        ci_d,ci_r,bad=m.bootstrap_effects(groups,np.array([1.]),np.random.default_rng(4),100)
        self.assertEqual(bad,100)
        self.assertTrue(np.isnan(ci_r).all())


class ReconciliationTests(unittest.TestCase):
    def test_legacy_gate_and_event_units(self):
        checks=pd.read_csv(m.OUT/"legacy_reproduction.csv")
        self.assertEqual(len(checks),20)
        self.assertTrue(checks.passed.all())
        e=pd.read_csv(m.OUT/"event_scores.csv")
        for v,g in e.groupby("virus_coordinate_set"):
            self.assertEqual(len(g),31201)
            self.assertEqual(g.unique_td_id.nunique(),31201)
            self.assertTrue((g.near == (g.left_hit | g.right_hit)).all())

    def test_zero_is_not_missing_and_fractions_reaggregate(self):
        d=pd.read_csv(m.OUT/"donor_scores.csv",float_precision="round_trip")
        self.assertTrue(d[d.true_zero_fraction].fraction.eq(0).all())
        self.assertTrue(d[d.total_DUP.isna()].fraction.isna().all())
        self.assertTrue((~d[d.total_DUP.isna()].included_corrected).all())
        e=pd.read_csv(m.OUT/"event_scores.csv")
        for (v,donor),g in e.groupby(["virus_coordinate_set","donor"]):
            row=d[(d.virus_coordinate_set==v)&(d.donor==donor)].iloc[0]
            self.assertEqual(row.near_DUP,int(g.near.sum()))
            self.assertEqual(row.total_DUP,len(g))
            self.assertAlmostEqual(row.fraction,g.near.mean(),places=11)

    def test_fraction_csv_roundtrip_retains_full_precision(self):
        d=pd.read_csv(m.OUT/"donor_scores.csv",float_precision="round_trip")
        valid=d[d.total_DUP.notna()]
        np.testing.assert_array_equal(valid.fraction.to_numpy(),(valid.near_DUP/valid.total_DUP).to_numpy())

    def test_brute_force_spatial_check(self):
        co=pd.read_csv(m.OUT/"coordinate_eligibility.csv")
        co=co[co.joint_eligible]
        e=pd.read_csv(m.OUT/"event_scores.csv").sample(80,random_state=44)
        for r in e.itertuples():
            ref=co[(co.virus_class==r.virus_coordinate_set)&(co.chrom==r.chrom)].pos0.to_numpy()
            manual=any(min(abs(int(p)-r.left_pos0),abs(int(p)-r.right_pos0))<=100000 for p in ref)
            self.assertEqual(manual,r.near)


if __name__=="__main__":
    unittest.main(verbosity=2)
