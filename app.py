import streamlit as st
import re
import pandas as pd

# ==========================================
# 1. 网页全局设置
# ==========================================
st.set_page_config(page_title="Fc 突变深度解码雷达 V17", page_icon="🛡️", layout="wide")

st.title("🛡️ 工业级 Fc 工程化突变解码雷达 (V17 智能排雷版)")
st.info("💡 终极指纹图谱：集成全球头部药企专利突变库。新增【同种异型 (Allotype) 鉴定】与【反向排雷 (Missing Mutation Alerts)】机制，预警 ADA 免疫原性风险与 CMC 纯化缺陷。")

# ==========================================
# 2. 核心知识库：全境 Fc 亚型、突变与同种异型图谱
# ==========================================
ISOTYPE_MOTIFS = {
    "IgG1": r"CPPCPAPELLG|CPPCPAPEAAG|CPPCPAPELAG",
    "IgG2": r"CVECPPCPAPPV",
    "IgG4 (野生型)": r"CPSCPAPEFLG",
    "IgG4 (S228P 稳定型)": r"CPPCPAPEFLG|CPPCPAPEAAG",
    "IgG3": r"CPRCPAPELLG",
}

# 同种异型鉴定字典 (Allotype)
ALLOTYPE_DB = {
    r"PPSR([D])E([L])T": ("G1m1", "D356/L358", "高加索人群常见; 强免疫原性标记"),
    r"PPSR([E])E([M])T": ("G1m3 / nG1m1", "E356/M358", "亚洲人群常见; G1m3 骨架标志"),
    r"VEP([K])SC": ("G1m17", "K214", "常与 G1m1 连锁 (CH1区)"),
    r"VEP([R])SC": ("nG1m17 (G1m3)", "R214", "常与 G1m3 连锁 (CH1区)"),
}

MUTATION_DB = {
    # --- 模块一：效应子功能沉默 ---
    r"PAPE([A])([A])GGPS": ("LALA", "L234A, L235A", "通用", "【毒性沉默】消除与 FcγR 结合，大幅降低 ADCC/CDC 效应"),
    r"PAPE([F])([E])GGPS": ("FEA", "L234F, L235E", "通用", "【毒性沉默】降低效应子功能 (常伴随D265A)"),
    r"PAPE([A])([A])G[A]P": ("PAA (IgG4)", "F234A, L235A", "通用", "【毒性沉默】IgG4专属的彻底静默突变"),
    r"VVV([S])VSHED": ("D265S", "D265S", "通用", "【毒性沉默】深度消除效应子结合，破坏结合界面"),
    r"VVV([A])VSHED": ("D265A", "D265A", "通用", "【毒性沉默】深度消除效应子结合，破坏结合界面"),
    r"NKAL([G])APIE": ("P329G", "P329G", "Roche", "【毒性沉默】破坏 FcγR 结合口袋，构成 LALA-PG 组合"),
    r"NNAKTKPREE[Q]Y[A]ST|NNAKTKPREE[Q]Y[Q]ST|NNAKTKPREE[Q]Y[G]ST": ("Aglycosylation", "N297A/Q/G", "通用", "【去糖基化】消除 N-糖基化位点，丧失效应子功能"),

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
    results = {"isotype": "未知 (未检测到标准骨架)", "mutations": [], "allotypes": [], "has_fc": False}
    
    if "LSPGK" in seq or "CPPCP" in seq or "VFLFPPKPK" in seq:
        results["has_fc"] = True
    else:
        return results

    # 检测亚型
    for iso, motif in ISOTYPE_MOTIFS.items():
        if re.search(motif, seq):
            results["isotype"] = iso
            break
            
    # 检测同种异型 (Allotypes)
    for motif, (allo_name, eu_num, desc) in ALLOTYPE_DB.items():
        if re.search(motif, seq):
            results["allotypes"].append({
                "名称": allo_name, "位置": eu_num, "描述": desc
            })

    # 检测特定专利突变
    for motif, (mut_name, eu_num, platform, func) in MUTATION_DB.items():
        if re.search(motif, seq):
            results["mutations"].append({
                "突变简称": mut_name, "EU 编号": eu_num, "技术溯源": platform, "生物学功能": func
            })
            
    return results

def parse_fasta(text):
    sequences = {}
    if ">" not in text:
        sequences["未命名输入序列_1"] = text.upper()
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
        deduction_reports = {} # 用于存储反向排雷所需的情报
        
        with st.spinner("正在遍历全球制药巨头突变专利库与同种异型字典..."):
            for name, seq in seq_dict.items():
                res = analyze_fc(seq)
                current_muts = [m["突变简称"] for m in res["mutations"]]
                current_allos = [a["名称"] for a in res["allotypes"]]
                
                # 保存情报，供智能战略模块使用
                deduction_reports[name] = {
                    "isotype": res["isotype"],
                    "muts": current_muts,
                    "allos": current_allos
                }
                
                if not res["has_fc"]:
                    report_data.append({"序列名称": name, "Fc 亚型骨架": "❌ 未检测到 Fc 区", "同种异型鉴定": "-", "特定突变识别": "-", "EU 编号": "-", "平台/溯源": "-", "临床/CMC 解析": "-"})
                    continue
                
                allo_str = " | ".join(current_allos) if current_allos else "未知或非典型"
                
                if not res["mutations"]:
                    report_data.append({"序列名称": name, "Fc 亚型骨架": res["isotype"], "同种异型鉴定": allo_str, "特定突变识别": "野生型 (WT)", "EU 编号": "-", "平台/溯源": "Natural", "临床/CMC 解析": "天然效应子功能与正常半衰期"})
                else:
                    for mut in res["mutations"]:
                        report_data.append({"序列名称": name, "Fc 亚型骨架": res["isotype"], "同种异型鉴定": allo_str, "特定突变识别": mut["突变简称"], "EU 编号": mut["EU 编号"], "平台/溯源": mut["技术溯源"], "临床/CMC 解析": mut["生物学功能"]})

        if report_data:
            df = pd.DataFrame(report_data)
            
            def highlight_rows(row):
                mut_str = str(row['特定突变识别'])
                if any(x in mut_str for x in ["Knob", "Hole", "EW", "RVT", "Azymetric", "Charge Steer"]): return ['background-color: #e3f2fd'] * len(row)
                if any(x in mut_str for x in ["LALA", "P329G", "D265", "FEA", "PAA"]): return ['background-color: #fff3e0'] * len(row)
                if any(x in mut_str for x in ["GA-SD", "AL-IE", "HexaBody"]): return ['background-color: #ffebee'] * len(row)
                if any(x in mut_str for x in ["Protein A", "S228P"]): return ['background-color: #e8f5e9'] * len(row)
                if any(x in mut_str for x in ["YTE", "LS", "IHH", "N434A"]): return ['background-color: #f3e5f5'] * len(row)
                return [''] * len(row)

            st.markdown("### 📊 Fc 全境工程化分析报告")
            st.dataframe(df.style.apply(highlight_rows, axis=1), use_container_width=True)
            
            # --- 智能战略级推演 & 反向排雷 (按序列隔离处理) ---
            st.markdown("### 💡 独立分子战略推演与反向排雷预警")
            
            for seq_name, data in deduction_reports.items():
                if data["isotype"] == "未知 (未检测到标准骨架)": continue
                
                muts = data["muts"]
                allos = data["allos"]
                iso = data["isotype"]
                
                with st.expander(f"📌 情报解密: {seq_name}", expanded=True):
                    
                    # ==========================================
                    # 💥 反向排雷预警 (Missing Mutation Alerts)
                    # ==========================================
                    has_warning = False
                    
                    # 1. IgG4 缺陷预警
                    if iso == "IgG4 (野生型)":
                        st.error("🚨 **反向排雷 [CMC风险]：缺失 S228P 稳定突变！** 检测到天然 IgG4 骨架。该分子在体内极易发生半分子交换 (Fab-arm exchange) 导致药物失效，强烈建议引入 S228P 或改成 IgG1 骨架。")
                        has_warning = True
                        
                    # 2. 双抗纯化缺陷预警
                    is_bispecific = any(x in "".join(muts) for x in ["Knob", "Hole", "EW", "Azymetric", "Charge Steer"])
                    if is_bispecific and "Protein A 破坏" not in muts:
                        st.warning("⚠️ **反向排雷 [下游工艺风险]：缺失不对称纯化突变！** 检测到双抗组装支架，但未见 H435R/Y436F。若无其他特殊纯化标签，同源二聚体杂质将极难通过常规 Protein A 层析去除。")
                        has_warning = True
                        
                    # 3. Allotype 嵌合冲突预警 (ADA 风险)
                    allo_str_joined = " ".join(allos)
                    if "G1m1" in allo_str_joined and "nG1m17" in allo_str_joined:
                        st.error("🚨 **反向排雷 [免疫原性风险]：同种异型冲突！** 同时检测到 G1m1 (D356/L358) 与 G1m3 特征 (R214)。这是一种非天然的人造嵌合体，可能具有极高的抗药抗体 (ADA) 激发风险。")
                        has_warning = True
                        
                    # 4. 突变负荷过高预警
                    if len(set(muts)) >= 4:
                        st.warning(f"⚠️ **反向排雷 [结构稳定性]：Fc 突变负荷过高 ({len(set(muts))}种)。** 叠加过多工程化突变极易导致 Fc 域微观折叠异常和新抗原表位暴露。")
                        has_warning = True

                    if not has_warning:
                        st.success("✅ **排雷扫描通过**：未见明显的 IgG4 缺失、纯化隐患或严重的同种异型冲突。")

                    st.markdown("---")
                    
                    # ==========================================
                    # 🎯 正向战略推演 (Strategic Deduction)
                    # ==========================================
                    st.markdown("##### 核心技术路线研判：")
                    # 1. 异源二聚化平台推演
                    if any("Azymetric" in m for m in muts): st.markdown("- 🎯 **双抗支架：Zymeworks (Azymetric) 不对称平台。**")
                    elif any("Charge Steer" in m for m in muts): st.markdown("- 🎯 **双抗支架：静电反转驱动 (Chugai/Xencor等)。**")
                    elif any("EW" in m for m in muts) or any("RVT" in m for m in muts): st.markdown("- 🎯 **双抗支架：高阶 EW-RVT 异源二聚化平台。**")
                    elif any("Knob" in m for m in muts) or any("Hole" in m for m in muts): st.markdown("- 🎯 **双抗支架：经典 Knob-in-Hole (Genentech)。**")
                        
                    # 2. 效应子功能推演
                    if "LALA" in muts and "P329G" in muts: st.markdown("- 🛡️ **杀伤机制：LALA-PG 终极效应子沉默。** (大概率为规避 CRS 的 T细胞接合器)")
                    elif "HexaBody" in muts: st.markdown("- ⚔️ **杀伤机制：补体风暴激发器 (Genmab HexaBody)。**")
                    elif any("GA-SD" in m for m in muts) or any("AL-IE" in m for m in muts): st.markdown("- ⚔️ **杀伤机制：超级 ADCC 增强。**")
                        
                    # 3. 药代动力学推演
                    if "YTE" in muts or "LS" in muts: st.markdown("- ⏱️ **PK 设计：超长效修饰。** (月度/季度给药设计)")
                    elif any("IHH" in m for m in muts): st.markdown("- ☢️ **PK 设计：极速体内清除。** (通常用于核药/ADC)")
                    
                    if not muts:
                        st.markdown("- 🧬 **常规抗体**：未检测到特殊的工程化修饰意图。")
    else:
        st.error("请输入序列！")
