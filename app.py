"""
DNAWorks Gradio Application
提供图形化界面调用Fortran版DNAWorks
"""

import gradio as gr
import subprocess
import os
import re
import tempfile
import pandas as pd
from pathlib import Path

# --- 配置区域 ---
DNAWORKS_PATH = "./dnaworks"
DEFAULT_CODON = "E. coli"
INP_DIR = "./inpfiles"
RESULT_DIR = "./results"

# 确保目录存在
Path(INP_DIR).mkdir(exist_ok=True)
Path(RESULT_DIR).mkdir(exist_ok=True)

# --- 预设密码子表 ---
# SerAGC会被错误识别为Asn，需要避免
# ProCCA仅会被唯一一个密码子识别，需要避免，ProCCU,ProCCC会造成+1位核糖体移码，应在重组基因中彻底避免
CODON_TABLES = {
    "E. coli": """CODOn
Gly     GGG     40359.00      11.39      0.16
Gly     GGA     34894.00       9.85      0.13
Gly     GGT     89915.00      25.37      0.35
Gly     GGC     94608.00      26.70      0.36
Glu     GAG     66665.00      18.81      0.33
Glu     GAA    137748.00      38.87      0.67
Asp     GAT    116164.00      32.78      0.63
Asp     GAC     67865.00      19.15      0.37
Val     GTG     85263.00      24.06      0.34
Val     GTA     41283.00      11.65      0.17
Val     GTT     70627.00      19.93      0.29
Val     GTC     50417.00      14.23      0.20
Ala     GCG    104293.00      29.43      0.32
Ala     GCA     75329.00      21.26      0.23
Ala     GCT     60787.00      17.15      0.19
Ala     GCC     85138.00      24.03      0.26
Arg     AGG      7966.00       2.25      0.04
Arg     AGA     13784.00       3.89      0.07
Ser     AGT     35966.00      10.15      0.16
Ser     AGC     53286.00      15.04      0.24
Lys     AAG     45133.00      12.74      0.26
Lys     AAA    125351.00      35.37      0.74
Asn     AAT     75086.00      21.19      0.50
Asn     AAC     75334.00      21.26      0.50
Met     ATG     92952.00      26.23      1.00
Ile     ATA     25982.00       7.33      0.12
Ile     ATT    105218.00      29.69      0.49
Ile     ATC     83118.00      23.46      0.39
Thr     ACG     48560.00      13.70      0.25
Thr     ACA     34483.00       9.73      0.17
Thr     ACT     37430.00      10.56      0.19
Thr     ACC     77023.00      21.74      0.39
Trp     TGG     48949.00      13.81      1.00
End     TGA      3616.00       1.02      0.31
Cys     TGT     18601.00       5.25      0.46
Cys     TGC     21434.00       6.05      0.54
End     TAG       978.00       0.28      0.08
End     TAA      7024.00       1.98      0.60
Tyr     TAT     62750.00      17.71      0.59
Tyr     TAC     43034.00      12.14      0.41
Leu     TTG     45581.00      12.86      0.13
Leu     TTA     51320.00      14.48      0.14
Phe     TTT     78743.00      22.22      0.58
Phe     TTC     56591.00      15.97      0.42
Ser     TCG     29993.00       8.46      0.13
Ser     TCA     32814.00       9.26      0.15
Ser     TCT     37586.00      10.61      0.17
Ser     TCC     32586.00       9.20      0.15
Arg     CGG     21391.00       6.04      0.11
Arg     CGA     13645.00       3.85      0.07
Arg     CGT     70009.00      19.76      0.36
Arg     CGC     68569.00      19.35      0.35
Gln     CAG    100346.00      28.32      0.66
Gln     CAA     51275.00      14.47      0.34
His     CAT     44633.00      12.60      0.58
His     CAC     32678.00       9.22      0.42
Leu     CTG    168885.00      47.66      0.47
Leu     CTA     15275.00       4.31      0.04
Leu     CTT     42704.00      12.05      0.12
Leu     CTC     35873.00      10.12      0.10
Pro     CCG     72450.00      20.44      0.49
Pro     CCA     30515.00       8.61      0.21
Pro     CCT     26805.00       7.56      0.18
Pro     CCC     19008.00       5.36      0.13
//""",
    "P. pastoris": """CODOn
Gly     GGG       468.00       5.76      0.00
Gly     GGA      1550.00      19.06      0.00
Gly     GGT      2075.00      25.52      0.00
Gly     GGC       655.00       8.06      0.00
Glu     GAG      2360.00      29.03      0.00
Glu     GAA      3043.00      37.43      0.00
Asp     GAT      2899.00      35.66      0.00
Asp     GAC      2103.00      25.87      0.00
Val     GTG       998.00      12.28      0.00
Val     GTA       804.00       9.89      0.00
Val     GTT      2188.00      26.91      0.00
Val     GTC      1210.00      14.88      0.00
Ala     GCG       314.00       3.86      0.00
Ala     GCA      1228.00      15.10      0.00
Ala     GCT      2351.00      28.92      0.00
Ala     GCC      1348.00      16.58      0.00
Arg     AGG       539.00       6.63      0.00
Arg     AGA      1634.00      20.10      0.00
Ser     AGT      1020.00      12.55      0.00
Ser     AGC       621.00       7.64      0.00
Lys     AAG      2748.00      33.80      0.00
Lys     AAA      2433.00      29.93      0.00
Asn     AAT      2038.00      25.07      0.00
Asn     AAC      2168.00      26.67      0.00
Met     ATG      1517.00      18.66      0.00
Ile     ATA       906.00      11.14      0.00
Ile     ATT      2532.00      31.14      0.00
Ile     ATC      1580.00      19.43      0.00
Thr     ACG       491.00       6.04      0.00
Thr     ACA      1118.00      13.75      0.00
Thr     ACT      1820.00      22.39      0.00
Thr     ACC      1175.00      14.45      0.00
Trp     TGG       834.00      10.26      0.00
End     TGA        27.00       0.33      0.00
Cys     TGT       626.00       7.70      0.00
Cys     TGC       356.00       4.38      0.00
End     TAG        40.00       0.49      0.00
End     TAA        69.00       0.85      0.00
Tyr     TAT      1300.00      15.99      0.00
Tyr     TAC      1473.00      18.12      0.00
Leu     TTG      2562.00      31.51      0.00
Leu     TTA      1265.00      15.56      0.00
Phe     TTT      1963.00      24.14      0.00
Phe     TTC      1675.00      20.60      0.00
Ser     TCG       598.00       7.36      0.00
Ser     TCA      1234.00      15.18      0.00
Ser     TCT      1983.00      24.39      0.00
Ser     TCC      1344.00      16.53      0.00
Arg     CGG       158.00       1.94      0.00
Arg     CGA       340.00       4.18      0.00
Arg     CGT       564.00       6.94      0.00
Arg     CGC       175.00       2.15      0.00
Gln     CAG      1323.00      16.27      0.00
Gln     CAA      2069.00      25.45      0.00
His     CAT       960.00      11.81      0.00
His     CAC       737.00       9.07      0.00
Leu     CTG      1215.00      14.94      0.00
Leu     CTA       873.00      10.74      0.00
Leu     CTT      1289.00      15.85      0.00
Leu     CTC       620.00       7.63      0.00
Pro     CCG       320.00       3.94      0.00
Pro     CCA      1540.00      18.94      0.00
Pro     CCT      1282.00      15.77      0.00
Pro     CCC       553.00       6.80      0.00
//""",
    "B. subtilis": """CODOn
Gly     GGG      9094.00     11.15      0.00
Gly     GGA     17743.00     21.76      0.00
Gly     GGT     10566.00     12.96      0.00
Gly     GGC     18967.00     23.26      0.00
Glu     GAG     18429.00     22.60      0.00
Glu     GAA     39217.00     48.09      0.00
Asp     GAT     27108.00     33.24      0.00
Asp     GAC     15485.00     18.99      0.00
Val     GTG     14105.00     17.30      0.00
Val     GTA     10638.00     13.05      0.00
Val     GTT     15193.00     18.63      0.00
Val     GTC     14067.00     17.25      0.00
Ala     GCG     16176.00     19.84      0.00
Ala     GCA     17205.00     21.10      0.00
Ala     GCT     15150.00     18.58      0.00
Ala     GCC     13433.00     16.47      0.00
Arg     AGG      3341.00      4.10      0.00
Arg     AGA      8561.00     10.50      0.00
Ser     AGT      5577.00      6.84      0.00
Ser     AGC     11766.00     14.43      0.00
Lys     AAG     16997.00     20.84      0.00
Lys     AAA     39449.00     48.38      0.00
Asn     AAT     18702.00     22.93      0.00
Asn     AAC     14522.00     17.81      0.00
Met     ATG     21424.00     26.27      0.00
Ile     ATA      7960.00      9.76      0.00
Ile     ATT     29487.00     36.16      0.00
Ile     ATC     22167.00     27.18      0.00
Thr     ACG     12110.00     14.85      0.00
Thr     ACA     17642.00     21.63      0.00
Thr     ACT      7082.00      8.68      0.00
Thr     ACC      7342.00      9.00      0.00
Trp     TGG      8765.00     10.75      0.00
End     TGA       621.00      0.76      0.00
Cys     TGT      2939.00      3.60      0.00
Cys     TGC      3472.00      4.26      0.00
End     TAG       381.00      0.47      0.00
End     TAA      1562.00      1.92      0.00
Tyr     TAT     18967.00     23.26      0.00
Tyr     TAC     10290.00     12.62      0.00
Leu     TTG     12914.00     15.84      0.00
Leu     TTA     16167.00     19.83      0.00
Phe     TTT     24450.00     29.98      0.00
Phe     TTC     11677.00     14.32      0.00
Ser     TCG      5266.00      6.46      0.00
Ser     TCA     11874.00     14.56      0.00
Ser     TCT     10320.00     12.66      0.00
Ser     TCC      6766.00      8.30      0.00
Arg     CGG      5664.00      6.95      0.00
Arg     CGA      3537.00      4.34      0.00
Arg     CGT      5903.00      7.24      0.00
Arg     CGC      6720.00      8.24      0.00
Gln     CAG     15089.00     18.50      0.00
Gln     CAA     16620.00     20.38      0.00
His     CAT     12832.00     15.74      0.00
His     CAC      6141.00      7.53      0.00
Leu     CTG     18757.00     23.00      0.00
Leu     CTA      3975.00      4.87      0.00
Leu     CTT     17772.00     21.79      0.00
Leu     CTC      8697.00     10.67      0.00
Pro     CCG     13316.00     16.33      0.00
Pro     CCA      5796.00      7.11      0.00
Pro     CCT      8640.00     10.60      0.00
Pro     CCC      2850.00      3.50      0.00
//"""
}


def validate_sequence(raw_input: str, mode: str) -> tuple:
    """
    验证并格式化输入序列
    Returns: (formatted_seq, clean_seq, error_msg)
    """
    # 清洗：移除非字母字符并转大写
    clean_seq = re.sub(r'[^a-zA-Z]', '', raw_input).upper()

    if not clean_seq:
        return "", "", "错误：输入序列不能为空或不包含有效字母！"

    if mode == 'PROTein':
        valid_chars = set("ACDEFGHIKLMNPQRSTVWY")
        label = "蛋白质"
        # 启发式拦截：核酸误贴检查
        nt_count = sum(1 for c in clean_seq if c in "ATCG")
        if len(clean_seq) > 10 and (nt_count / len(clean_seq)) > 0.9:
            return "", "", "风险拦截：检测到当前极像是【核酸序列】，但您选择了【蛋白质模式】！"
    else:
        valid_chars = set("ATCG")
        label = "核酸"

    invalid_found = [c for c in clean_seq if c not in valid_chars]
    if invalid_found:
        illegal_str = " ".join(set(invalid_found))
        return "", "", f"错误：{label}模式下包含非法字符: {illegal_str}"

    # DNAWorks 格式：60字符一行，带行号
    lines = []
    for i in range(0, len(clean_seq), 60):
        lines.append(f"      {i+1:<3} {clean_seq[i:i+60]}")
    formatted_str = "\n".join(lines)

    return formatted_str, clean_seq, ""


def generate_inp_content(
    title: str,
    mode: str,
    clean_seq: str,
    formatted_seq: str,
    melting: int,
    len_low: int,
    len_high: int,
    frequency: int,
    codon_table: str,
    log_path: str
) -> str:
    """生成 DNAWorks 输入文件内容"""

    codon_section = CODON_TABLES.get(codon_table, codon_table) if mode == 'PROTein' else ""

    content = f'''title "{title}"
melting low {melting}
length low {len_low} high {len_high}
frequency threshold {frequency}
logfile "{log_path}"

PATTern
   BsaI GGTCTC
   BspQI GCTCTTC
   SacI GAGCTC
   SalI GTCGAC
   BglII AGATCT
   Esp3I CGTCTC
//

{codon_section}

{mode}
{formatted_seq}
//
'''
    return content


def run_dnaworks(
    title: str,
    mode: str,
    sequence: str,
    melting: int,
    len_low: int,
    len_high: int,
    frequency: int,
    codon_table: str,
    custom_codon_table: str = ""
) -> tuple:
    """
    主运行函数
    Returns: (status_msg, result_text, primer_df_html)
    """
    # 验证序列
    formatted_seq, clean_seq, error = validate_sequence(sequence, mode)
    if error:
        return f"❌ {error}", "", ""

    # 生成文件名
    name = title.strip().replace(' ', '_')
    suffix = "_pro" if mode == 'PROTein' else "_nuc"
    inp_path = os.path.join(INP_DIR, f"{name}{suffix}.inp")
    log_path = os.path.join(RESULT_DIR, f"{name}{suffix}.txt")

    # 生成输入文件
    # 如果选择了"自定义"，则使用用户输入的自定义密码子表
    effective_codon_table = custom_codon_table if codon_table == "自定义" else codon_table
    inp_content = generate_inp_content(
        title, mode, clean_seq, formatted_seq,
        melting, len_low, len_high, frequency, effective_codon_table, log_path
    )

    try:
        with open(inp_path, 'w') as f:
            f.write(inp_content)
    except Exception as e:
        return f"❌ 写入输入文件失败: {str(e)}", "", ""

    # 检查可执行文件
    if not os.path.exists(DNAWORKS_PATH):
        return f"❌ 错误：在路径 '{DNAWORKS_PATH}' 找不到可执行文件。请先编译 Fortran 代码。", "", ""

    # 运行 DNAWorks
    try:
        result = subprocess.run(
            [DNAWORKS_PATH, inp_path],
            capture_output=True,
            text=True,
            timeout=300
        )
        stdout = result.stdout
        stderr = result.stderr

        if result.returncode != 0:
            return f"❌ 程序异常退出，退出码: {result.returncode}\n{stderr}", "", ""

    except subprocess.TimeoutExpired:
        return "❌ 运行超时（5分钟）", "", ""
    except Exception as e:
        return f"❌ 系统错误: {str(e)}", "", ""

    # 读取结果
    if not os.path.exists(log_path):
        return "⚠️ DNAWorks  运行完成，但未找到结果文件", stdout, ""

    try:
        with open(log_path, 'r') as f:
            log_content = f.read()
    except Exception as e:
        return f"❌ 读取结果文件失败: {str(e)}", "", ""

    # 解析引物
    primer_html = extract_primers_html(log_content)

    return f"✅ 运行完成！\n\n--- DNAWorks 输出 ---\n{stdout}", log_content, primer_html


def extract_primers_html(log_path: str, trial_num: int = 1) -> str:
    """从日志文件中提取引物表格（HTML格式）"""
    try:
        with open(log_path, 'r') as f:
            content = f.read()
    except:
        return "<p>无法读取结果文件</p>"

    # 定位 Trial 区域
    pattern = rf"PARAMETERS FOR TRIAL\s+{trial_num}\b(.*?)(?:PARAMETERS FOR TRIAL|$)"
    trial_section = re.search(pattern, content, re.DOTALL | re.IGNORECASE)

    if not trial_section:
        return f"<p>未找到 Trial {trial_num} 的数据</p>"

    section_text = trial_section.group(1)

    # 匹配引物行: 名称,序列,长度
    primer_pattern = r"(\w+),([ATCG]+),(\d+)nt"
    matches = re.findall(primer_pattern, section_text)

    if not matches:
        return "<p>该 Trial 区域内未识别到引物列表</p>"

    data = [{"引物名称": p[0], "序列 (5'→3')": p[1], "长度 (nt)": p[2]} for p in matches]
    df = pd.DataFrame(data)

    return df.to_html(index=False, escape=False, classes="primer-table")


# --- Gradio 界面 ---
css = """
.gradio-container {max-width: 1200px !important}
h1 {text-align: center; color: #2c3e50;}
.primer-table {font-size: 14px; border-collapse: collapse; width: 100%;}
.primer-table th {background-color: #3498db; color: white; padding: 8px; text-align: left;}
.primer-table td {padding: 8px; border-bottom: 1px solid #ddd;}
.primer-table tr:hover {background-color: #f5f5f5;}
"""

with gr.Blocks(title="DNAWorks 寡核苷酸设计", css=css) as demo:
    gr.Markdown("# 🧬 DNAWorks 寡核苷酸设计工具")
    gr.Markdown("*基于 PCR 的基因合成引物自动化设计*")

    with gr.Row():
        with gr.Column(scale=1):
            gr.Markdown("### 📝 基本参数")

            title = gr.Textbox(
                label="任务名称",
                value="INS_HUMAN",
                info="用于命名输出文件"
            )

            mode = gr.Radio(
                label="输入类型",
                choices=[("蛋白质序列 (Protein)", "PROTein"), ("核酸序列 (Nucleotide)", "NUCLeotide")],
                value="PROTein",
                info="选择输入的是蛋白序列还是核酸序列"
            )

            sequence = gr.Textbox(
                label="序列",
                placeholder="粘贴蛋白序列或核酸序列...",
                value="FVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
                info="支持直接粘贴，包含数字和空格的序列会自动清洗",
                lines=6
            )

        with gr.Column(scale=1):
            gr.Markdown("### ⚙️ 高级参数")

            melting = gr.Slider(
                label="第二轮 Tm (℃)",
                minimum=50,
                maximum=80,
                value=62,
                step=1,
                info="第一轮 Tm 会自动设置为比此值低 4~5℃"
            )

            len_low = gr.Slider(
                label="引物最短长度 (nt)",
                minimum=45,
                maximum=59,
                value=50,
                step=1
            )

            len_high = gr.Slider(
                label="引物最长长度 (nt)",
                minimum=45,
                maximum=59,
                value=57,
                step=1
            )

            codon_table = gr.Dropdown(
                label="密码子优化物种",
                choices=list(CODON_TABLES.keys()) + ["自定义"],
                value="E. coli",
                info="仅对蛋白质模式有效"
            )

            custom_codon_table = gr.Textbox(
                label="自定义密码子表",
                placeholder="请按照以下格式输入密码子表...\nCODOn\nGly     GGG     40359.00      11.39      0.16\n...",
                lines=10,
                visible=False,
                info="当上方选择「自定义」时，此处填写自定义密码子表"
            )

            frequency = gr.Slider(
                label="最低密码子频率 (%)",
                minimum=5,
                maximum=30,
                value=20,
                step=1,
                info="仅使用频率高于此值的密码子"
            )

    run_btn = gr.Button("🚀 生成并运行 DNAWorks", variant="primary", size="lg")

    gr.Markdown("---")

    with gr.Row():
        with gr.Column(scale=1):
            gr.Markdown("### 📊 运行状态")
            status_output = gr.Textbox(label="状态", lines=3, interactive=False)

        with gr.Column(scale=1):
            gr.Markdown("### 📋 结果预览 (Trial 1)")
            result_output = gr.Textbox(label="完整日志", lines=15, interactive=False)

    gr.Markdown("---")

    with gr.Row():
        gr.Markdown("### 🧪 引物列表 (Trial 1)")
        primer_output = gr.HTML(label="引物表格")

    # 事件绑定
    def update_custom_input(selection):
        return gr.update(visible=(selection == "自定义"))

    codon_table.change(
        fn=update_custom_input,
        inputs=[codon_table],
        outputs=[custom_codon_table]
    )

    run_btn.click(
        fn=run_dnaworks,
        inputs=[title, mode, sequence, melting, len_low, len_high, frequency, codon_table, custom_codon_table],
        outputs=[status_output, result_output, primer_output]
    )

    # 实时更新第一轮Tm提示
    melting.change(
        fn=lambda x: f"ℹ️ 第一轮 Tm 将设置为: {x - 4}℃",
        inputs=[melting],
        outputs=[]
    )

# 启动
if __name__ == "__main__":
    demo.launch()
