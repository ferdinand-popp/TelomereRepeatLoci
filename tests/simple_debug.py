# test/run_cli.py
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from debug_wrapper import make_debug_args, run_direct_pipeline  # noqa: E402


def main():
    # Make sure we run from the project root (where `src/` and `data/` live)
    project_root = PROJECT_ROOT
    args = make_debug_args(
        tumor_bam=str(
            project_root
            / "data/HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram"
        ),
        tel_tumor_bam=str(
            project_root
            / "data/HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage/tumor_TelomerCnt_HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage/HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage_filtered_intratelomeric.bam"
        ),
        control_bam=str(
            project_root
            / "data/HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram"
        ),
        tel_control_bam=str(
            project_root
            / "data/HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage/tumor_TelomerCnt_HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage/HG00152.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage_filtered_intratelomeric.bam"
        ),
        # downloaded with curl -L -C - http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -o GRCh38_full_analysis_set_plus_decoy_hla.fa
        reference_fasta=str(
            project_root / "data/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        ),
    )
    print("Running direct pipeline wrapper", file=sys.stderr)
    run_direct_pipeline(args)


if __name__ == "__main__":
    main()
