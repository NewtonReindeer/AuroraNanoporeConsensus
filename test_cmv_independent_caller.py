import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).parent
spec = importlib.util.spec_from_file_location('cmv_independent_caller', ROOT / 'cmv_independent_caller.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

class IndependentTests(unittest.TestCase):
    def test_optout_does_not_touch_inputs(self):
        obj = {'cmv_variant_calls': [{'foo': 'bar'}]}
        with patch.dict(os.environ, {'AURORA_CMV_INDEPENDENT_CALLER': '0'}):
            self.assertIs(module.optional_run('s', '/missing.fastq', '/missing', '/missing', obj), obj)

    def test_reference_coordinate_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp)/'r.fa'
            fa.write_text('>NC_006273.2\nACCGTT\n')
            self.assertEqual(module.validate_reference(fa)['length'], 6)
            vcf = Path(tmp)/'x.vcf'
            vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nNC_006273.2\t2\t.\tC\tT\t60\tPASS\t.\n')
            self.assertEqual(module.validate_vcf(vcf, fa), 1)
            vcf.write_text(vcf.read_text().replace('\t2\t.\tC\tT', '\t2\t.\tA\tT'))
            with self.assertRaises(module.ReviewError):
                module.validate_vcf(vcf, fa)

    def test_disabled_or_missing_config_is_review_only(self):
        with patch.dict(os.environ, {'AURORA_CMV_INDEPENDENT_CALLER': '1', 'AURORA_CMV_MERLIN_FASTA': '', 'AURORA_CMV_CLAIR3_MODEL': ''}):
            result = module.optional_run('s', 'missing.fastq', '.', '.', {'cmv_variant_calls': [{'a':1}]})
            self.assertTrue(result['cmv_independent_status'].startswith('REVIEW_FAILED:'))
            self.assertEqual(result['cmv_variant_calls'], [{'a':1}])

    def test_reference_contig_lock(self):
        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp)/'r.fa'
            fa.write_text('>other\nACGT\n')
            with self.assertRaises(module.ReviewError):
                module.validate_reference(fa)

if __name__ == '__main__':
    unittest.main()
