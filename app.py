import streamlit as st
import re
import pandas as pd

# ==========================================
# 1. 网页全局设置
# ==========================================
st.set_page_config(page_title="Fc 突变深度解码雷达", page_icon="🛡️", layout="wide")

st.title("🛡️ 工业级 Fc 工程化突变解码雷达 (全平台统治版)")
st.info("💡 终极指纹图谱：已集成 Genentech (KIH)、Zymeworks (Azymetric)、Xencor (XmAb)、Chugai (ART-Ig)、Genmab (HexaBody) 等全球领先药企的核心专利突变群。")

# ==========================================
# 2. 核心知识库：全境 Fc 亚型与突变指纹图谱
# ==========================================
ISOTYPE_MOTIFS = {
    "IgG1": r"CPPCPAPELLG|CPPCPAPEAAG|CPPCPAPELAG",
    "IgG2": r"CVECPPCPAPPV",
    "IgG4 (野生型)": r"CPSCPAPEFLG",
    "IgG4 (S228P 稳定型)": r"CPPCPAPEFLG|CPPCPAPEAAG",
    "IgG3": r"CPRCPAPELLG",
}

MUTATION_DB = {
    # --- 模块一：效应子功能沉默 ---
    r"PAPE([A])([A])GGPS": ("LALA", "L234A, L235A", "通用", "【毒性沉默】消除与 FcγR 结合，大幅降低 ADCC/CDC 效应"),
    r"PAPE([F])([E])GGPS": ("FEA", "L234F, L235E", "通用", "【毒性沉默】降低效应子功能 (常伴随D265A)"),
    r"PAPE([A])([A])G[A]P": ("PAA (IgG4)", "F234A, L235A", "通用", "【毒性沉默】IgG4专属的彻底静默突变"),
    r"VVV([S])VSHED": ("D265S", "D265S", "通用", "【毒性沉默】深度消除效应子结合，破坏结合界面"),
    r"VVV([A])VSHED": ("D265A", "D265A", "通用", "【毒性沉默】深度消除效应子结合，破坏结合界面"),
    r"NKAL([G])APIE": ("P329G", "P329G", "Roche", "【毒性沉默】破坏 FcγR 结合口袋，构成 LALA-PG 组合"),
    r"NNAKTKPREE[Q]Y[A]ST|NNAKTKPREE[Q]Y[Q]ST|NNAKTKPREE[Q]Y[G]ST": ("Aglycosylation", "N297A/Q/G", "通用", "【去糖基化】消除 N-糖基化位点，导致完全丧失效应子功能"),

    # --- 模块二：效应子功能增强 ---
    r"PAPELL([A])GP([D])VFL": ("GA-SD", "G236A, S239D", "Xencor", "【功能增强】显著增强对 FcγRIIa 和 FcγRIIIa 的亲和力"),
    r"NKAL([L])P([E])EKT": ("AL-IE", "A330L, I332E", "Xencor", "【功能增强】配合 S239D 构成 GASDALIE，实现 ADCC 极限增强"),
    r"TQKSLSLS([G])PGK": ("HexaBody", "E430G", "Genmab", "【六聚体化】在细胞表面促进形成六聚体，引爆极限 CDC"),
    r"REEM([R])NQVS": ("HexaBody", "E345R", "Genmab", "【六聚体化】同 E430G，促进寡聚化以增强 CDC"),

    # --- 模块三：经典与高阶异源二聚体化 ---
    r"NQVSL([W])CLVK": ("Knob (凸起)", "T366W", "Genentech", "【双抗组装】经典 KIH 凸起端，强制异源二聚化"),
    r"NQVSL([S])C([A])VK": ("Hole 核心", "T366S, L368A", "Genentech", "【双抗组装】经典 KIH 凹陷端，容纳 Knob"),
    r"DGSFFL([V])SKLT": ("Hole 尾部", "Y407V", "Genentech", "【双抗组装】经典 KIH 凹陷端，扩大疏水口袋"),
    
    r"REEMT([E])NQVSL": ("EW 平台 (E)", "K392E", "EW-RVT", "【双抗组装】电荷反转，与 Chain B 的 R 配对"),
    r"FFLYS([W])LTVD": ("EW 平台 (W)", "K409W", "EW-RVT", "【双抗组装】引入疏水凸起，插入 Chain B 凹陷"),
    r"PREP([R])VYTL": ("RVT 平台 (R)", "Q347R", "EW-RVT", "【双抗组装】引入正电荷，与 Chain A 的 E 配对"),
    r"PPVL([V])SDGSF([T])LYS": ("RVT 平台 (VT)", "D399V, F405T", "EW-RVT", "【双抗组装】构建疏水容纳凹陷"),

    r"VYTLPP([V])([Y])REEMT": ("Azymetric Chain A (VY)", "T350V, L351Y", "Zymeworks", "【双抗组装】Azymetric 不对称支架 A 链前段"),
    r"SDGS([A])L([V])SKLT": ("Azymetric Chain A (AV)", "F405A, Y407V", "Zymeworks", "【双抗组装】Azymetric 不对称支架 A 链尾段"),
    r"VYTLPP([V])SREEM": ("Azymetric Chain B (V)", "T350V", "Zymeworks", "【双抗组装】Azymetric 不对称支架 B 链"),
    r"NQVSL([L])CLVK": ("Azymetric Chain B (L)", "T366L", "Zymeworks", "【双抗组装】Azymetric 不对称支架 B 链核心"),
    r"PENNY([L])T([W])PPVL": ("Azymetric Chain B (LW)", "K392L, T394W", "Zymeworks", "【双抗组装】Azymetric 不对称支架 B 链空间匹配"),

    r"REEMT([K])NQVS": ("Charge Steer Chain A (K)", "E356K", "Chugai/XmAb", "【双抗组装】负转正电荷，利用静电排斥同源二聚"),
    r"PPVL([K])SDGS": ("Charge Steer Chain A (K)", "D399K", "Chugai/XmAb", "【双抗组装】负转正电荷，深度静电互补"),
    r"PENNY([D])TTPP": ("Charge Steer Chain B (D)", "K392D", "Chugai/XmAb", "【双抗组装】正转负电荷，配合 Chain A 的 K 突变"),
    r"FFLYS([D])LTV": ("Charge Steer Chain B (D)", "K409D", "Chugai/XmAb", "【双抗组装】正转负电荷，配合 Chain A 的 K 突变"),

    # --- 模块四：药代动力学优化 ---
    r"DTL([Y])I([T])R([E])PE": ("YTE", "M252Y, S254T, T256E", "AstraZeneca", "【PK优化】增强酸性下对 FcRn 的亲和力，延长半衰期"),
    r"HEALHNH([L])TQ([S])": ("LS", "M428L, N434S", "Xencor", "【PK优化】增强 FcRn 亲和力，当前最主流的长效修饰方案"),
    r"HEALHNH([A])TQK": ("N434A", "N434A", "通用", "【PK优化】适度增强 FcRn 结合，延长半衰期"),
    r"PKDTL([A])ISRT": ("IHH 减半衰期 (I)", "I253A", "通用", "【快速清除】彻底破坏 FcRn 结合，极速缩短半衰期"),
    r"HQDWL([A])GKEY": ("IHH 减半衰期 (H1)", "H310A", "通用", "【快速清除】协同 I253A，阻止抗体被再循环"),

    # --- 模块五：CMC 工艺与纯化适配 ---
    r"HEALHN([R])([F])TQK": ("Protein A 破坏", "H435R, Y436F", "通用", "【CMC纯化】破坏与 Protein A 的结合，利用梯度洗脱纯化双抗"),
    r"CPPCPAPEFLG": ("S228P", "S228P", "通用", "【结构稳定】IgG4 专属，强制形成链间二硫键，防止 Fab 臂交换")
}

# ==========================================
# 3. 核心扫描引擎
# ==========================================
def analyze_fc(sequence):
    seq = sequence.upper().replace(" ", "").replace("\n", "")
    results = {"isotype": "未知 (未检测到标准骨架)", "mutations": [], "has_fc": False}
    
    if "LSPGK" in seq or "CPPCP" in seq or "VFLFPPKPK" in seq:
        results["has_fc"] = True
    else:
        return results

    for iso, motif in ISOTYPE_MOTIFS.items():
        if re.search(motif, seq):
            results["isotype"] = iso
            break
            
    for motif, (mut_name, eu_num, platform, func) in MUTATION_DB.items():
        if re.search(motif, seq):
            results["mutations"].append({
                "突变简称": mut_name,
                "EU 编号": eu_num,
                "技术溯源": platform,
                "生物学功能": func
            })
    return results

def parse_fasta(text):
    sequences = {}
    if ">" not in text:
        sequences["未命名输入序列"] = text.upper()
        return sequences
    for part in text.split(">"):
        if not part.strip(): continue
        lines = part.strip().split("\n")
        name = lines[0].strip()
        seq = "".join(lines[1:]).strip().upper()
        if name and seq: sequences[name] = seq
    return sequences

# ==========================================
# 4. 交互界面
# ==========================================
raw_input = st.text_area("📥 粘贴抗体全长链或 Fc 段序列 (支持多条 FASTA):", height=200)

if st.button("🔍 启动全境 Fc 深度解码", type="primary"):
    if raw_input:
        seq_dict = parse_fasta(raw_input)
        report_data = []
        deduction_reports = {} # 独立存储每条序列的推演情报
        
        with st.spinner("正在遍历全球制药巨头突变专利库..."):
            for name, seq in seq_dict.items():
                res = analyze_fc(seq)
                current_seq_muts = []
                
                if not res["has_fc"]:
                    report_data.append({"序列名称": name, "Fc 亚型骨架": "❌ 未检测到 Fc 区", "检测到的特定突变": "-", "EU 编号": "-", "技术溯源": "-", "临床/CMC 功能解析": "-"})
                    continue
                if not res["mutations"]:
                    report_data.append({"序列名称": name, "Fc 亚型骨架": res["isotype"], "检测到的特定突变": "野生型 (WT)", "EU 编号": "-", "技术溯源": "Natural", "临床/CMC 功能解析": "天然效应子功能与正常半衰期"})
                else:
                    for mut in res["mutations"]:
                        current_seq_muts.append(mut["突变简称"])
                        report_data.append({"序列名称": name, "Fc 亚型骨架": res["isotype"], "检测到的特定突变": mut["突变简称"], "EU 编号": mut["EU 编号"], "技术溯源": mut["技术溯源"], "临床/CMC 功能解析": mut["生物学功能"]})
                
                # 记录该序列的独有突变组合，供下游独立推演
                deduction_reports[name] = current_seq_muts

        if report_data:
            df = pd.DataFrame(report_data)
            
            def highlight_rows(row):
                mut_str = str(row['检测到的特定突变'])
                if any(x in mut_str for x in ["Knob", "Hole", "EW", "RVT", "Azymetric", "Charge Steer"]): return ['background-color: #e3f2fd'] * len(row)
                if any(x in mut_str for x in ["LALA", "P329G", "D265", "FEA", "PAA"]): return ['background-color: #fff3e0'] * len(row)
                if any(x in mut_str for x in ["GA-SD", "AL-IE", "HexaBody"]): return ['background-color: #ffebee'] * len(row)
                if any(x in mut_str for x in ["Protein A", "S228P"]): return ['background-color: #e8f5e9'] * len(row)
                if any(x in mut_str for x in ["YTE", "LS", "IHH", "N434A"]): return ['background-color: #f3e5f5'] * len(row)
                return [''] * len(row)

            st.markdown("### 📊 Fc 全境工程化分析报告")
            st.dataframe(df.style.apply(highlight_rows, axis=1), use_container_width=True)
            
            # --- 智能战略级推演 (按序列隔离处理) ---
            st.markdown("### 💡 独立分子战略意图推演")
            
            for seq_name, muts in deduction_reports.items():
                if not muts: continue
                
                with st.expander(f"📌 情报解密: {seq_name}", expanded=True):
                    # 1. 异源二聚化平台推演
                    if any("Azymetric" in m for m in muts):
                        st.success("🎯 **支架判定：Zymeworks (Azymetric) 不对称平台。** 高技术壁垒的疏水/空间联合突变，用于极高纯度的双抗组装。")
                    elif any("Charge Steer" in m for m in muts):
                        st.success("🎯 **支架判定：静电反转驱动双抗。** 依靠电荷相吸的技术，通常来自 Chugai (ART-Ig) 或 Xencor。")
                    elif any("EW" in m for m in muts) or any("RVT" in m for m in muts):
                        st.success("🎯 **支架判定：高阶 EW-RVT 异源二聚化平台。** 成功检测到 EWRVT 配对网络特征。")
                    elif any("Knob" in m for m in muts) or any("Hole" in m for m in muts):
                        st.success("🎯 **支架判定：经典 Knob-in-Hole (Genentech) 双抗。** 最正统的凸起与凹陷空间位阻设计。")
                        
                    # 2. 效应子功能推演
                    if "LALA" in muts and "P329G" in muts:
                        st.warning("⚠️ **杀伤判定：LALA-PG 终极效应子沉默。** 业界最严苛的去毒性组合，极大概率为规避 CRS 的 T细胞接合器 (TCE)。")
                    elif "HexaBody" in muts:
                        st.error("🔥 **杀伤判定：补体风暴激发器 (Genmab HexaBody)。** 诱导抗体形成六聚体，产生毁灭性的 CDC 杀伤效力。")
                    elif any("GA-SD" in m for m in muts) or any("AL-IE" in m for m in muts):
                        st.error("🔥 **杀伤判定：超级 ADCC 增强 (GASDALIE等)。** 显著提升巨噬细胞/NK细胞招募，剑指高强度实体瘤杀伤。")
                        
                    # 3. 药代动力学推演
                    if "YTE" in muts or "LS" in muts:
                        st.info("⏱️ **PK 判定：超长效修饰。** 药物被设计为每月甚至每季度一次给药的长效制剂。")
                    elif any("IHH" in m for m in muts):
                        st.info("☢️ **PK 判定：极速体内清除。** 故意破坏 FcRn 结合，这通常是 ADC 毒素载体、放射性核素偶联药物 (RDC) 的标志！")
    else:
        st.error("请输入序列！")
