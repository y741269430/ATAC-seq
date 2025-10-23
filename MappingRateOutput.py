import os
import re
import sys
import pandas as pd

def parse_mapping_output(file_path):
    """解析单个 HISAT2/Bowtie2 输出文件，提取关键统计指标"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    results = {'filename': os.path.basename(file_path)}
    
    # 总 reads 数
    total_reads = re.search(r'(\d+) reads; of these:', content)
    if total_reads:
        results['total_reads'] = int(total_reads.group(1))
    
    # 完全不对齐的 reads（concordantly 0 次）
    concordantly_0 = re.search(r'(\d+) \(([\d.]+)%\) aligned concordantly 0 times', content)
    if concordantly_0:
        results['concordantly_0_times'] = int(concordantly_0.group(1))
        results['concordantly_0_times_%'] = float(concordantly_0.group(2))
    
    # 完全对齐 1 次的 reads
    concordantly_1 = re.search(r'(\d+) \(([\d.]+)%\) aligned concordantly exactly 1 time', content)
    if concordantly_1:
        results['concordantly_exactly_1_time'] = int(concordantly_1.group(1))
        results['concordantly_exactly_1_time_%'] = float(concordantly_1.group(2))
    
    # 对齐超过 1 次的 reads
    concordantly_gt1 = re.search(r'(\d+) \(([\d.]+)%\) aligned concordantly >1 times', content)
    if concordantly_gt1:
        results['concordantly_>1_times'] = int(concordantly_gt1.group(1))
        results['concordantly_>1_times_%'] = float(concordantly_gt1.group(2))
    
    # 不一致对齐的 reads（discordantly 1 次）
    discordantly_1 = re.search(r'(\d+) \(([\d.]+)%\) aligned discordantly 1 time', content)
    if discordantly_1:
        results['discordantly_1_time'] = int(discordantly_1.group(1))
        results['discordantly_1_time_%'] = float(discordantly_1.group(2))
    
    # 总体对齐率
    overall_alignment = re.search(r'([\d.]+)% overall alignment rate', content)
    if overall_alignment:
        results['overall_alignment_rate_%'] = float(overall_alignment.group(1))
    
    return results

def process_mapinfo_directory(directory, output_basename):
    """批量处理目录中的 HISAT2 输出文件，并输出 CSV 和 TXT 文件"""
    data = []
    
    for filename in os.listdir(directory):
        if filename.endswith('_map.txt'):
            file_path = os.path.join(directory, filename)
            try:
                data.append(parse_mapping_output(file_path))
            except Exception as e:
                print(f"警告: 解析文件 {file_path} 时出错: {e}")
                continue
    
    if not data:
        print(f"错误: 在目录 {directory} 中没有找到任何以 '_map.txt' 结尾的文件")
        sys.exit(1)
    
    # 构建 DataFrame
    df = pd.DataFrame(data)
    
    # 定义期望的列顺序（可以根据你的实际输出调整）
    columns_order = [
        'filename', 'total_reads',
        'concordantly_0_times', 'concordantly_0_times_%',
        'concordantly_exactly_1_time', 'concordantly_exactly_1_time_%',
        'concordantly_>1_times', 'concordantly_>1_times_%',
        'discordantly_1_time', 'discordantly_1_time_%',
        'overall_alignment_rate_%'
    ]
    
    # 只保留 DataFrame 中实际存在的列
    columns_order = [col for col in columns_order if col in df.columns]
    df = df[columns_order]
    
    # 输出文件名
    csv_filename = output_basename + '.csv'
    txt_filename = output_basename + '.txt'
    
    # 保存为 CSV
    df.to_csv(csv_filename, index=False)
    print(f"✅ CSV 结果已保存到: {csv_filename}")
    
    # 保存为 TXT（制表符分隔）
    df.to_csv(txt_filename, index=False, sep='\t')
    print(f"✅ TXT（制表符分隔）结果已保存到: {txt_filename}")
    
    return df

if __name__ == "__main__":
    # 检查命令行参数
    if len(sys.argv) != 3:
        print("❌ 用法: python MappingRate.py 输入目录 输出文件名（不含扩展名）")
        print("✅ 示例: python MappingRate.py mapinfo alignment_results")
        print(" 脚本将生成: alignment_results.csv 和 alignment_results.txt")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_basename = sys.argv[2]  # 不要带 .csv 或 .txt
    
    if not os.path.isdir(input_dir):
        print(f"❌ 错误: 输入目录 '{input_dir}' 不存在")
        sys.exit(1)
    
    try:
        df = process_mapinfo_directory(input_dir, output_basename)
        print("\n 整合后的数据预览:")
        print(df.head())
    except Exception as e:
        print(f"❌ 处理过程中发生错误: {e}")
        sys.exit(1)