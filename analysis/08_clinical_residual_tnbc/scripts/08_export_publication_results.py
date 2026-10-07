"""Export supplied existing clinical result tables without fitting models or scoring genes."""
import argparse
from pathlib import Path
from clinical_io import rows,save

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-dir',required=True);p.add_argument('--output-dir',required=True);a=p.parse_args()
    for name in ['SourceData_Seo_clinical_patient_scores.tsv','SourceData_Seo_clinical_model_results.tsv','SourceData_Seo_program146_mapping.tsv','SourceData_Seo_incremental_model.tsv']:
        save(Path(a.output_dir)/name,rows(Path(a.input_dir)/name))
if __name__=='__main__':main()
