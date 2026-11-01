#!/usr/bin/env python3
"""
add_exogenous_genes.py

将外源基因序列添加到参考基因组 FASTA 文件和 GTF 文件中。
输出文件在原文件名基础上加 .add_exo 后缀。
"""

import argparse
import shutil
from pathlib import Path


def add_genes(reference_fasta, reference_gtf, exo_fasta, source="external"):
    # 生成输出文件名
    out_fasta = str(Path(reference_fasta).with_suffix(".add_exo.fasta"))
    out_gtf = str(Path(reference_gtf).with_suffix(".add_exo.gtf"))

    # 1. 拷贝原始文件
    shutil.copy(reference_fasta, out_fasta)
    shutil.copy(reference_gtf, out_gtf)
    print(f"已复制原始 FASTA -> {out_fasta}，GTF -> {out_gtf}")

    # 2. 读取外源基因序列
    genes = {}
    with open(exo_fasta, "r") as f:
        seq_name = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_name:
                    genes[seq_name] = "".join(seq_lines)
                seq_name = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_name:
            genes[seq_name] = "".join(seq_lines)

    # 3. 追加外源基因到 FASTA（整条序列一行）
    with open(out_fasta, "a") as f:
        for gene, seq in genes.items():
            f.write(f">{gene}\n")
            f.write(seq + "\n")

    # 4. 生成 GTF 条目并追加
    with open(out_gtf, "a") as f:
        for gene, seq in genes.items():
            length = len(seq)
            transcript = f"{gene}-001"
            f.write(f'{gene}\t{source}\tgene\t1\t{length}\t.\t+\t.\tgene_id "{gene}"; gene_name "{gene}";\n')
            f.write(
                f'{gene}\t{source}\ttranscript\t1\t{length}\t.\t+\t.\tgene_id "{gene}"; transcript_id "{transcript}";\n'
            )
            f.write(
                f'{gene}\t{source}\texon\t1\t{length}\t.\t+\t.\tgene_id "{gene}"; transcript_id "{transcript}"; exon_number "1";\n'
            )

    print(f"成功将 {len(genes)} 个外源基因添加到 {out_fasta} 和 {out_gtf}")


def main():
    parser = argparse.ArgumentParser(description="将外源基因添加到参考 FASTA 和 GTF 中")
    parser.add_argument("--ref_fasta", required=True, help="原始参考基因组 FASTA")
    parser.add_argument("--ref_gtf", required=True, help="原始 GTF 注释文件")
    parser.add_argument("--exo_fasta", required=True, help="外源基因 FASTA 文件")
    parser.add_argument("--source", default="external", help="GTF source 字段 (default: external)")
    args = parser.parse_args()

    add_genes(args.ref_fasta, args.ref_gtf, args.exo_fasta, args.source)


if __name__ == "__main__":
    main()

