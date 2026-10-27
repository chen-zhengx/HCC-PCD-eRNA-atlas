#!/bin/bash
# =========================================================
# 对应结果：表达预处理第 1 步——fastp 质控与剪切
# 用法：bash 01_fastp_batch.sh <原始 FASTQ 目录> <输出目录>   （FASTQ 命名为 <SRR>_1.fastq.gz / <SRR>_2.fastq.gz）
# =========================================================
# 设置目录
INPUT_DIR="${1:?usage: $0 INPUT_DIR OUTPUT_DIR}"
OUTPUT_DIR="${2:?usage: $0 INPUT_DIR OUTPUT_DIR}"
mkdir -p "$OUTPUT_DIR"
# 日志文件配置
PROCESS_LOG="$OUTPUT_DIR/processing_details.log"    # 处理进度日志
FASTP_REPORT_LOG="$OUTPUT_DIR/fastp_reports.log"    # 所有fastp报告汇总
TRUE_ERROR_LOG="$OUTPUT_DIR/real_failures.log"      # 仅记录真实错误
SUMMARY_LOG="$OUTPUT_DIR/summary.log"               # 处理汇总日志
# 初始化日志文件
echo "【处理进度日志】- $(date)" > "$PROCESS_LOG"
echo "【fastp报告汇总】- $(date)" > "$FASTP_REPORT_LOG"
echo "【真实错误记录】- $(date)" > "$TRUE_ERROR_LOG"
echo "【处理汇总】- $(date)" > "$SUMMARY_LOG"
# 终端显示（带颜色）
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
NC='\033[0m' # No Color
# 获取所有样本名称（通过查找所有_1.fastq.gz文件）
echo -e "${YELLOW}⏳ 扫描输入目录中的样本文件...${NC}"
SAMPLES=($(find "$INPUT_DIR" -name "*_1.fastq.gz" -exec basename {} _1.fastq.gz \; | sort))
if [ ${#SAMPLES[@]} -eq 0 ]; then
    error_msg="[ERROR] 在输入目录中未找到任何样本文件: $INPUT_DIR"
    echo -e "${RED}✗ $error_msg${NC}"
    echo "$error_msg" >> "$TRUE_ERROR_LOG"
    exit 1
fi
echo -e "${GREEN}找到 ${#SAMPLES[@]} 个样本: ${SAMPLES[*]}${NC}"
echo "找到 ${#SAMPLES[@]} 个样本: ${SAMPLES[*]}" >> "$PROCESS_LOG"
# 初始化计数器
total_count=${#SAMPLES[@]}
success_count=0
fail_count=0
# 遍历所有样本
for SAMPLE in "${SAMPLES[@]}"; do
    R1="$INPUT_DIR/${SAMPLE}_1.fastq.gz"
    R2="$INPUT_DIR/${SAMPLE}_2.fastq.gz"
    echo -e "\n${BLUE}▶${NC} 正在处理: ${BOLD}$SAMPLE${NC} ($((success_count + fail_count + 1))/$total_count)"
    # 记录处理进度
    echo "[$(date +'%F %T')] 开始处理: $SAMPLE" >> "$PROCESS_LOG"
    echo "输入文件: $R1" >> "$PROCESS_LOG"
    echo "         $R2" >> "$PROCESS_LOG"
    # 检查文件是否存在
    if [[ ! -f "$R1" ]]; then
        error_msg="[ERROR] 文件缺失: $SAMPLE 的 _1.fastq.gz 不存在"
        echo -e "${RED}✗ $error_msg${NC}"
        echo "$error_msg" >> "$TRUE_ERROR_LOG"
        echo "$error_msg" >> "$PROCESS_LOG"
        ((fail_count++))
        continue
    fi
    if [[ ! -f "$R2" ]]; then
        error_msg="[ERROR] 文件缺失: $SAMPLE 的 _2.fastq.gz 不存在"
        echo -e "${RED}✗ $error_msg${NC}"
        echo "$error_msg" >> "$TRUE_ERROR_LOG"
        echo "$error_msg" >> "$PROCESS_LOG"
        ((fail_count++))
        continue
    fi
    # 检查文件是否可读
    if [[ ! -r "$R1" || ! -r "$R2" ]]; then
        error_msg="[ERROR] 文件不可读: $SAMPLE 的fastq文件权限不足"
        echo -e "${RED}✗ $error_msg${NC}"
        echo "$error_msg" >> "$TRUE_ERROR_LOG"
        echo "$error_msg" >> "$PROCESS_LOG"
        ((fail_count++))
        continue
    fi
    # 运行fastp（将报告输出到汇总文件）
    echo -e "\n===== $SAMPLE fastp报告 =====" >> "$FASTP_REPORT_LOG"
    echo -e "${YELLOW}⏳ 运行fastp质控...${NC}"
    start_time=$(date +%s)
    if fastp -i "$R1" -I "$R2" \
       -o "$OUTPUT_DIR/${SAMPLE}_1.clean.fq.gz" \
       -O "$OUTPUT_DIR/${SAMPLE}_2.clean.fq.gz" \
       -h "$OUTPUT_DIR/${SAMPLE}.html" \
       -j "$OUTPUT_DIR/${SAMPLE}.json" \
       -q 20 -n 15 -l 36 -w 8 >> "$FASTP_REPORT_LOG" 2>> "$TRUE_ERROR_LOG"
    then
        end_time=$(date +%s)
        duration=$((end_time - start_time))
        # 成功处理
        echo -e "${GREEN}✓ 完成: $SAMPLE (耗时: ${duration}s)${NC}"
        echo "[$(date +'%F %T')] 成功完成: $SAMPLE (耗时: ${duration}s)" >> "$PROCESS_LOG"
        # 检查输出文件是否生成
        if [[ -f "$OUTPUT_DIR/${SAMPLE}_1.clean.fq.gz" && -f "$OUTPUT_DIR/${SAMPLE}_2.clean.fq.gz" ]]; then
            echo "输出文件: ${SAMPLE}_1.clean.fq.gz, ${SAMPLE}_2.clean.fq.gz" >> "$PROCESS_LOG"
            ((success_count++))
        else
            error_msg="[WARNING] fastp运行成功但输出文件可能未正确生成: $SAMPLE"
            echo -e "${YELLOW}⚠ $error_msg${NC}"
            echo "$error_msg" >> "$PROCESS_LOG"
            ((fail_count++))
        fi
    else
        end_time=$(date +%s)
        duration=$((end_time - start_time))
        # 处理失败
        error_msg="[ERROR] fastp处理失败: $SAMPLE (退出码: $? 耗时: ${duration}s)"
        echo -e "${RED}✗ $error_msg${NC}"
        echo "$error_msg" >> "$TRUE_ERROR_LOG"
        echo "$error_msg" >> "$PROCESS_LOG"
        # 清理可能生成的不完整文件
        rm -f "$OUTPUT_DIR/${SAMPLE}_1.clean.fq.gz" "$OUTPUT_DIR/${SAMPLE}_2.clean.fq.gz" \
               "$OUTPUT_DIR/${SAMPLE}.html" "$OUTPUT_DIR/${SAMPLE}.json"
        ((fail_count++))
    fi
    echo "----------------------------------" >> "$PROCESS_LOG"
done
# 最终汇总
echo -e "\n${BOLD}=== 处理结果汇总 ===${NC}"
echo -e "${GREEN}成功样本: $success_count${NC}"
if [[ $fail_count -gt 0 ]]; then
    echo -e "${RED}失败样本: $fail_count${NC}"
else
    echo -e "${GREEN}失败样本: $fail_count${NC}"
fi
echo -e "总共样本: $total_count"
# 记录到汇总日志
echo "=== 处理结果汇总 ===" >> "$SUMMARY_LOG"
echo "处理时间: $(date)" >> "$SUMMARY_LOG"
echo "输入目录: $INPUT_DIR" >> "$SUMMARY_LOG"
echo "输出目录: $OUTPUT_DIR" >> "$SUMMARY_LOG"
echo "总共样本: $total_count" >> "$SUMMARY_LOG"
echo "成功样本: $success_count" >> "$SUMMARY_LOG"
echo "失败样本: $fail_count" >> "$SUMMARY_LOG"
echo "成功率: $(echo "scale=2; $success_count * 100 / $total_count" | bc)%" >> "$SUMMARY_LOG"
echo "详细日志: $PROCESS_LOG" >> "$SUMMARY_LOG"
echo "fastp报告: $FASTP_REPORT_LOG" >> "$SUMMARY_LOG"
if [[ $fail_count -gt 0 ]]; then
    echo -e "${RED}错误日志: $TRUE_ERROR_LOG${NC}"
    echo "错误日志: $TRUE_ERROR_LOG" >> "$SUMMARY_LOG"
fi
# 显示生成的输出文件统计
echo -e "\n${BOLD}输出文件统计:${NC}"
echo "清理后的fastq文件: $(find "$OUTPUT_DIR" -name "*clean.fq.gz" | wc -l)"
echo "HTML报告文件: $(find "$OUTPUT_DIR" -name "*.html" | wc -l)"
echo "JSON报告文件: $(find "$OUTPUT_DIR" -name "*.json" | wc -l)"
# 显示处理完成的样本列表
if [ $success_count -gt 0 ]; then
    echo -e "\n${GREEN}成功处理的样本:${NC}"
    find "$OUTPUT_DIR" -name "*clean.fq.gz" | sed 's/.*\///' | sort | uniq | sed 's/_[12]\.clean\.fq\.gz//' | uniq
fi
if [ $fail_count -gt 0 ]; then
    echo -e "\n${RED}处理失败的样本:${NC}"
    grep "ERROR" "$TRUE_ERROR_LOG" | grep -o "SRR[0-9]*" | sort | uniq
fi
