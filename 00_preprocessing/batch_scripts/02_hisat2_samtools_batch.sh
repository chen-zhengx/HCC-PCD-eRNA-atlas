#!/bin/bash
# =========================================================
# 对应结果：表达预处理第 2 步——HISAT2 比对、samtools 排序与建索引、比对统计
# 用法：bash 02_hisat2_samtools_batch.sh <fastp 输出目录> <输出目录> <HISAT2 索引前缀>
# =========================================================
# ====== 用户配置区域 ======
input_dir="${1:?usage: $0 INPUT_DIR OUTPUT_DIR HISAT2_INDEX}"
output_dir="${2:?usage: $0 INPUT_DIR OUTPUT_DIR HISAT2_INDEX}"
index="${3:?usage: $0 INPUT_DIR OUTPUT_DIR HISAT2_INDEX}"
threads=8  # 使用线程数，可以根据服务器资源调整
# ====== 错误检查函数 ======
check_error() {
    local status=$1
    local step="$2"
    local sample="$3"
    local logfile="$4"
    if [ $status -ne 0 ]; then
        echo "❌ [ERROR] 样本 ${sample} 的步骤【${step}】失败！" | tee -a "$logfile"
        return 1
    fi
    return 0
}
# ====== 主程序开始 ======
mkdir -p "${output_dir}" || { echo "❌ 输出目录创建失败！"; exit 1; }
# 创建主日志文件
master_log="${output_dir}/master.log"
echo "🎯 开始批量处理样本 - $(date)" | tee "$master_log"
echo "输入目录: ${input_dir}" | tee -a "$master_log"
echo "输出目录: ${output_dir}" | tee -a "$master_log"
echo "索引路径: ${index}" | tee -a "$master_log"
echo "线程数: ${threads}" | tee -a "$master_log"
# 获取所有样本名称（通过查找所有_1.clean.fq.gz文件）
echo "🔍 扫描输入目录中的样本文件..." | tee -a "$master_log"
samples=($(find "$input_dir" -name "*_1.clean.fq.gz" -exec basename {} _1.clean.fq.gz \; | sort))
if [ ${#samples[@]} -eq 0 ]; then
    echo "❌ 在输入目录中未找到任何样本文件: $input_dir" | tee -a "$master_log"
    echo "查找的文件模式: *_1.clean.fq.gz" | tee -a "$master_log"
    exit 1
fi
total_samples=${#samples[@]}
echo "✅ 找到 ${total_samples} 个样本: ${samples[*]}" | tee -a "$master_log"
# 初始化计数器
success_count=0
fail_count=0
failed_samples=()
# 遍历所有样本
for ((i=0; i<${#samples[@]}; i++)); do
    sample="${samples[i]}"
    current_num=$((i + 1))
    echo -e "\n📋 处理样本 ${current_num}/${total_samples}: ${sample}" | tee -a "$master_log"
    # 为每个样本创建单独的日志文件
    logfile="${output_dir}/${sample}.log"
    echo "🎯 开始处理样本：${sample} - $(date)" > "$logfile"
    echo "输入目录: ${input_dir}" >> "$logfile"
    echo "输出目录: ${output_dir}" >> "$logfile"
    # 检查输入文件是否存在
    if [[ ! -f "${input_dir}/${sample}_1.clean.fq.gz" ]]; then
        error_msg="❌ 输入文件不存在: ${input_dir}/${sample}_1.clean.fq.gz"
        echo "$error_msg" | tee -a "$master_log" >> "$logfile"
        ((fail_count++))
        failed_samples+=("$sample")
        continue
    fi
    if [[ ! -f "${input_dir}/${sample}_2.clean.fq.gz" ]]; then
        error_msg="❌ 输入文件不存在: ${input_dir}/${sample}_2.clean.fq.gz"
        echo "$error_msg" | tee -a "$master_log" >> "$logfile"
        ((fail_count++))
        failed_samples+=("$sample")
        continue
    fi
    # 比对 - 使用管道直接生成排序的BAM文件
    echo "[1/3] 正在进行 HISAT2 比对并直接生成排序BAM..." | tee -a "$logfile"
    hisat2 -x "${index}" \
        -1 "${input_dir}/${sample}_1.clean.fq.gz" \
        -2 "${input_dir}/${sample}_2.clean.fq.gz" \
        -p ${threads} \
        --dta \
        --rna-strandness RF \
        --new-summary \
        --summary-file "${output_dir}/${sample}.hisat2.summary" \
        2>>"$logfile" \
        | samtools sort -@ ${threads} -o "${output_dir}/${sample}.sorted.bam" -
    if check_error $? "HISAT2比对" "$sample" "$logfile"; then
        # 创建BAM索引
        echo "[2/3] 创建BAM文件索引..." | tee -a "$logfile"
        samtools index -@ ${threads} "${output_dir}/${sample}.sorted.bam"
        if check_error $? "BAM索引创建" "$sample" "$logfile"; then
            # 生成统计信息
            echo "[3/3] 生成统计信息..." | tee -a "$logfile"
            {
                echo "===== HISAT2 比对报告 ====="
                cat "${output_dir}/${sample}.hisat2.summary"
                echo -e "\n===== BAM 文件统计 ====="
                samtools stats "${output_dir}/${sample}.sorted.bam" | \
                grep -E '^(SN|ID)' | \
                awk -F'\t' '$2 ~ /reads mapped:|reads properly paired:/ {print $2 $3}'
                echo -e "\n===== 比对率统计 ====="
                samtools flagstat "${output_dir}/${sample}.sorted.bam" | grep -E "mapped|properly"
            } > "${output_dir}/${sample}.full.report" 2>>"$logfile"
            if check_error $? "统计信息生成" "$sample" "$logfile"; then
                echo "✅ 样本 ${sample} 处理完成！" | tee -a "$master_log" >> "$logfile"
                echo "📊 统计文件: ${output_dir}/${sample}.full.report" >> "$logfile"
                echo "🔧 比对文件: ${output_dir}/${sample}.sorted.bam" >> "$logfile"
                ((success_count++))
            else
                ((fail_count++))
                failed_samples+=("$sample")
            fi
        else
            ((fail_count++))
            failed_samples+=("$sample")
        fi
    else
        ((fail_count++))
        failed_samples+=("$sample")
    fi
    echo "----------------------------------" >> "$logfile"
done
# 最终汇总报告
echo -e "\n\n📈 ===== 批量处理完成汇总 =====" | tee -a "$master_log"
echo "🕒 完成时间: $(date)" | tee -a "$master_log"
echo "✅ 成功样本数: ${success_count}/${total_samples}" | tee -a "$master_log"
echo "❌ 失败样本数: ${fail_count}/${total_samples}" | tee -a "$master_log"
if [ ${#failed_samples[@]} -gt 0 ]; then
    echo "📋 失败样本列表: ${failed_samples[*]}" | tee -a "$master_log"
fi
# 显示输出文件统计
echo -e "\n💾 输出文件统计:" | tee -a "$master_log"
echo "BAM文件: $(find "${output_dir}" -name "*.sorted.bam" | wc -l)" | tee -a "$master_log"
echo "BAM索引: $(find "${output_dir}" -name "*.bai" | wc -l)" | tee -a "$master_log"
echo "统计报告: $(find "${output_dir}" -name "*.full.report" | wc -l)" | tee -a "$master_log"
echo "HISAT2摘要: $(find "${output_dir}" -name "*.hisat2.summary" | wc -l)" | tee -a "$master_log"
# 显示磁盘使用情况
echo -e "\n📊 磁盘使用情况:" | tee -a "$master_log"
du -sh "${output_dir}" | tee -a "$master_log"
if [ $fail_count -eq 0 ]; then
    echo -e "\n🎉 所有样本处理成功完成！" | tee -a "$master_log"
else
    echo -e "\n⚠️  有 ${fail_count} 个样本处理失败，请检查日志文件。" | tee -a "$master_log"
fi
echo "📝 主日志文件: ${master_log}" | tee -a "$master_log"
