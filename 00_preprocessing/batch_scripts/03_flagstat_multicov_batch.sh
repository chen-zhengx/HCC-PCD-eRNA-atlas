#!/bin/bash
# =========================================================
# 对应结果：表达预处理第 3 步——samtools flagstat（RPM 分母）与 bedtools multicov（eRNA 区域计数）
# 输出：<flagstat 目录>/<SRR>_aligned_sorted.txt、<计数目录>/<SRR>_bd.txt；之后由 ../code/01_erna_rpm.R 合并为 count 与 RPM 矩阵
# 用法：bash 03_flagstat_multicov_batch.sh <BAM 目录> <eRNA bed> <flagstat 输出目录> <计数输出目录>
#       （BAM 命名为 <SRR>_sorted.bam，并有 .bai 索引）
# =========================================================
# 配置路径
pbam_dir="${1:?usage: $0 BAM_DIR BED_FILE FLAGSTAT_DIR COUNT_DIR}"
bed_file="${2:?usage: $0 BAM_DIR BED_FILE FLAGSTAT_DIR COUNT_DIR}"
wy_output="${3:?usage: $0 BAM_DIR BED_FILE FLAGSTAT_DIR COUNT_DIR}"
bd_output="${4:?usage: $0 BAM_DIR BED_FILE FLAGSTAT_DIR COUNT_DIR}"
# 创建日志目录（按时间戳命名）
log_dir="./logs/$(date '+%Y%m%d_%H%M%S')"
global_log="$log_dir/global.log"
mkdir -p "$wy_output" "$bd_output" "$log_dir"
# 初始化统计变量
total_files=$(ls "$pbam_dir"/*_sorted.bam | grep -v '.bai$' | wc -l)
processed=0
success_count=0
declare -a fail_records
# 获取有效BAM文件列表
mapfile -t bam_files < <(ls "$pbam_dir"/*_sorted.bam | grep -v '.bai$')
# 打印初始化信息
echo "===== 开始处理 $total_files 个BAM文件 =====" | tee -a "$global_log"
echo "日志目录：$log_dir" | tee -a "$global_log"
echo "当前时间：$(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$global_log"
# 主处理循环
for bam in "${bam_files[@]}"
do
    ((processed++))
    base=$(basename "$bam" _sorted.bam)
    file_log="$log_dir/${base}.log"
    index_file="$pbam_dir/${base}_sorted.bam.bai"
    echo -e "\n===== 处理进度 $processed/$total_files：$base =====" | tee -a "$global_log"
    echo "开始时间：$(date '+%H:%M:%S')" | tee -a "$file_log" "$global_log"
    # 检查索引文件
    if [ ! -f "$index_file" ]; then
        err_msg="错误：缺失索引文件 ${base}_sorted.bam.bai"
        echo "$err_msg" | tee -a "$file_log" "$global_log"
        fail_records+=("$base - 致命错误：$err_msg")
        continue
    fi
    # 步骤1：samtools flagstat
    echo "[1/2] 正在运行samtools flagstat..." | tee -a "$file_log"
    if samtools flagstat "$bam" > "$wy_output/${base}_aligned_sorted.txt" 2>> "$file_log"
    then
        echo "√ flagstat完成" | tee -a "$file_log"
    else
        err_msg="flagstat执行失败，退出码：$?"
        echo "$err_msg" | tee -a "$file_log" "$global_log"
        fail_records+=("$base - 步骤1失败：$err_msg")
        continue
    fi
    # 步骤2：bedtools multicov
    echo "[2/2] 正在运行bedtools multicov..." | tee -a "$file_log"
    if bedtools multicov -bams "$bam" -bed "$bed_file" > "$bd_output/${base}_bd.txt" 2>> "$file_log"
    then
        echo "√ multicov完成" | tee -a "$file_log"
        ((success_count++))
    else
        err_msg="multicov执行失败，退出码：$?"
        echo "$err_msg" | tee -a "$file_log" "$global_log"
        fail_records+=("$base - 步骤2失败：$err_msg")
    fi
    echo "结束时间：$(date '+%H:%M:%S')" | tee -a "$file_log" "$global_log"
done
# 生成统计报告
echo -e "\n===== 处理结果统计 =====" | tee -a "$global_log"
echo "成功文件数：$success_count/$total_files" | tee -a "$global_log"
if [ ${#fail_records[@]} -gt 0 ]; then
    echo -e "\n===== 错误详情 =====" | tee -a "$global_log"
    printf "%-15s %-s\n" "样本编号" "错误描述" | tee -a "$global_log"
    for record in "${fail_records[@]}"; do
        IFS=" - " read -r id msg <<< "$record"
        printf "%-15s %-s\n" "$id" "$msg" | tee -a "$global_log"
    done
fi
echo -e "\n详细日志查看：" | tee -a "$global_log"
echo "全局日志：$global_log" | tee -a "$global_log"
echo "个体日志：$log_dir/[SRR编号].log" | tee -a "$global_log"
