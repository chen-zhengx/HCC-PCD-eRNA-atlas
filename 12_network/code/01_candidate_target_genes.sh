#!/bin/bash
# =========================================================
# 对应结果：网络中风险 eRNA 的候选靶基因（eRNA 上下游 ±1 Mb 范围内的基因；Methods "Construction of eRNA networks"）
# 输入：NCBI RefSeq 人类基因组注释 GCF_000001405.40（GRCh38.p14）的 genomic.gtf（可从 NCBI Datasets 下载）
#       input/risk_eRNAs_by_PCD.csv（12 个风险 eRNA 的位置及其染色体的 RefSeq 编号）
# 输出：output/target_genes/<eRNA>/genes_<范围>bp.txt（基因名、起点、终点；±100 kb 至 ±1 Mb）与 summary_gene_counts.txt
# 运行：在模块目录下 bash code/01_candidate_target_genes.sh /path/to/genomic.gtf（awk 的三参数 match() 需要 GNU awk）
# =========================================================

gtf="$1"
if [ -z "$gtf" ] || [ ! -f "$gtf" ]; then
    echo "用法：bash code/01_candidate_target_genes.sh /path/to/genomic.gtf"
    exit 1
fi

tail -n +2 input/risk_eRNAs_by_PCD.csv | tr -d '\r' | cut -d, -f3,5,6,7 | sort -u | while IFS=, read -r erna start end chr; do

output_dir="output/target_genes/${erna//:/}"
mkdir -p "$output_dir"
echo "== ${erna} (${chr}:${start}-${end})"

# 查询范围 ±100kb, ±200kb, ±300kb, ±400kb, ±500kb, ±1000kb
ranges=(100000 200000 300000 400000 500000 1000000)

# 统计文件
summary_file="${output_dir}/summary_gene_counts.txt"
echo -e "Range(bp)\tGene_Count" > "$summary_file"

for range in "${ranges[@]}"; do
    outfile="${output_dir}/genes_${range}bp.txt"
    echo "查询范围 ±$((range/1000))kb 的基因..."
    awk -v chr="$chr" -v start=$start -v end=$end -v range=$range \
    '$1==chr && $3=="gene" && $4 <= (end+range) && $5 >= (start-range) {
        if (match($0, /gene "([^"]+)"/, arr)) {
            print arr[1], $4, $5;
        }
    }' "$gtf" | sort -u > "$outfile"
    count=$(wc -l < "$outfile")
    echo -e "±$((range/1000))kb\t$count" >> "$summary_file"
done

echo "统计文件已生成：$summary_file"

done
