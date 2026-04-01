import os
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed


def parse_args():
    parser = argparse.ArgumentParser(description="Run GCTA heritability analysis by group.")
    parser.add_argument("--group", required=True, help="Group name, e.g. CA or CT")
    parser.add_argument("--gcta", default="gcta64", help="Path to GCTA executable")
    parser.add_argument("--bfile-prefix", required=True, help="PLINK binary prefix")
    parser.add_argument("--phe-file", required=True, help="Phenotype file")
    parser.add_argument("--qcovar-file", required=True, help="Quantitative covariate file")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--trait-start", type=int, default=1, help="Start phenotype index")
    parser.add_argument("--trait-end", type=int, default=4491, help="End phenotype index")
    parser.add_argument("--grm-threads", type=int, default=30, help="Threads for GRM construction")
    parser.add_argument("--reml-threads", type=int, default=10, help="Threads for each REML run")
    parser.add_argument("--max-workers", type=int, default=10, help="Parallel worker number")
    return parser.parse_args()


def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)


def build_grm(args, grm_prefix):
    grm_bin = f"{grm_prefix}.grm.bin"
    grm_nbin = f"{grm_prefix}.grm.N.bin"
    grm_id = f"{grm_prefix}.grm.id"

    if os.path.isfile(grm_bin) and os.path.isfile(grm_nbin) and os.path.isfile(grm_id):
        print("[INFO] GRM files already exist, skipping GRM construction")
        return

    cmd_grm = (
        f"{args.gcta} --bfile {args.bfile_prefix} "
        f"--autosome-num 35 --autosome "
        f"--make-grm --make-grm-alg 1 "
        f"--thread-num {args.grm_threads} "
        f"--out {grm_prefix}"
    )
    print(f"[INFO] Building GRM: {cmd_grm}")
    run_command(cmd_grm)


def run_h2(idx, trait, args, grm_prefix):
    trait_prefix = os.path.join(args.output_dir, f"h2_{idx}_{trait}")
    hsq_file = f"{trait_prefix}.hsq"

    if not os.path.isfile(hsq_file):
        cmd_h2 = (
            f"{args.gcta} --reml --grm {grm_prefix} "
            f"--reml-pred-rand "
            f"--pheno {args.phe_file} "
            f"--qcovar {args.qcovar_file} "
            f"--mpheno {idx} "
            f"--thread-num {args.reml_threads} "
            f"--out {trait_prefix}"
        )
        print(f"[RUN] {trait}: {cmd_h2}")
        run_command(cmd_h2)

    if os.path.isfile(hsq_file):
        with open(hsq_file) as fo:
            for line in fo:
                if "V(G)/Vp" in line:
                    parts = line.strip().split("\t")
                    h2_val = round(float(parts[1]), 4)
                    h2_se = round(float(parts[2]), 4)
                    h2_std = f"{round(float(parts[1]), 2)}±{round(float(parts[2]), 2)}"
                    return f"{trait}\t{h2_val}\t{h2_se}\t{h2_std}\n"

    return None


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    phe_list = [f"meta{i}" for i in range(args.trait_start, args.trait_end + 1)]
    grm_prefix = os.path.join(args.output_dir, f"{args.group}_Vanraden_grm1")
    summary_file = os.path.join(args.output_dir, "h2_summary.txt")

    build_grm(args, grm_prefix)

    results = []
    print(f"[INFO] Starting parallel REML analysis with {args.max_workers} workers")

    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        future_to_trait = {
            executor.submit(run_h2, idx, trait, args, grm_prefix): trait
            for idx, trait in enumerate(phe_list, start=args.trait_start)
        }

        for future in as_completed(future_to_trait):
            trait = future_to_trait[future]
            try:
                res = future.result()
                if res:
                    results.append(res)
                    print(f"[DONE] {trait}")
            except Exception as e:
                print(f"[ERROR] {trait}: {e}")

    with open(summary_file, "w") as sum_w:
        sum_w.writelines(results)

    print(f"[INFO] All traits finished. Summary saved to: {summary_file}")


if __name__ == "__main__":
    main()