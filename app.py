import streamlit as st
import pandas as pd
import re
from Bio import Align
import py3Dmol
import streamlit.components.v1 as components

# ==========================================
# 1. 网页全局设置
# ==========================================
st.set_page_config(page_title="Fc 突变深度解码雷达 V19.2", page_icon="🛡️", layout="wide")

st.title("🛡️ 工业级 Fc 工程化突变解码雷达 (V19.2 终极稳定版)")
st.info("💡 终极形态：【Biopython 精准空间对齐】+【满血版彩色战略情报大屏】+【原生 HTML 注入固化的 3D 映射实验室】。")

# ==========================================
# 2. 核心知识库：野生型标尺与空间坐标字典
# ==========================================
WT_START_EU = 214
WT_REF_SEQ = "KREPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"

ISOTYPE_MOTIFS = {
    "IgG1": r"CPPCPAPELLG|CPPCPAPEAAG|CPPCPAPELAG",
    "IgG2": r"CVECPPCPAPPV",
    "IgG4 (野生型)": r"CPSCPAPEFLG",
    "IgG4 (S228P 稳定型)": r"CPPCPAPEFLG|CPPCPAPEAAG",
    "IgG3": r"CPRCPAPELLG",
}

ALLOTYPE_DB = {
    "G1m1": ({"356": 'D', "358": 'L'}, "D356/L358", "高加索人群常见; 强免疫原性标记"),
    "G1m3 / nG1m1": ({"356": 'E', "358": 'M'}, "E356/M358", "亚洲人群常见; G1m3 骨架标志"),
    "G1m17": ({"214": 'K'}, "K214", "常与 G1m1 连锁 (CH1区)"),
    "nG1m17 (G1m3)": ({"214": 'R'}, "R214", "常与 G1m3 连锁 (CH1区)"),
}

MUTATION_DB = {
    "LALA": ({"234": 'A', "235": 'A'}, "L234A, L235A", "通用", "【毒性沉默】消除与 FcγR 结合，大幅降低 ADCC/CDC 效应"),
    "PAA (IgG4)": ({"234": 'F', "235": 'A'}, "F234A, L235A", "通用", "【毒性沉默】IgG4专属的彻底静默突变"),
    "FEA": ({"234": 'F', "235": 'E'}, "L234F, L235E", "通用", "【毒性沉默】降低效应子功能"),
    "D265S": ({"265": 'S'}, "D265S", "通用", "【毒性沉默】深度消除效应子结合，破坏结合界面"),
    "D265A": ({"265": 'A'}, "D265A", "通用", "【毒性沉默】深度消除效应子结合，破坏结合界面"),
    "P329G": ({"329": 'G'}, "P329G", "Roche", "【毒性沉默】破坏 FcγR 结合口袋，构成 LALA-PG 组合"),
    "Aglycosylation": ({"297": ['A', 'Q', 'G']}, "N297A/Q/G", "通用", "【去糖基化】消除 N-糖基化位点，丧失效应子功能"),

    "GA-SD": ({"236": 'A', "239": 'D'}, "G236A, S239D", "Xencor", "【功能增强】显著增强对 FcγRIIa 和 FcγRIIIa 的亲和力"),
    "AL-IE": ({"330": 'L', "332": 'E'}, "A330L, I332E", "Xencor", "【功能增强】配合 S239D 构成 GASDALIE，实现 ADCC 极限增强"),
    "HexaBody (E430G)": ({"430": 'G'}, "E430G", "Genmab", "【六聚体化】在细胞表面促进形成六聚体，引爆极限 CDC"),
    "HexaBody (E345R)": ({"345": 'R'}, "E345R", "Genmab", "【六聚体化】同 E430G，促进寡聚化以增强 CDC"),

    "Knob (凸起)": ({"366": 'W'}, "T366W", "Genentech", "【双抗组装】经典 KIH 凸起端，强制异源二聚化"),
    "Hole (核心)": ({"366": 'S', "368": 'A'}, "T366S, L368A", "Genentech", "【双抗组装】经典 KIH 凹陷端，容纳 Knob"),
    "Hole (尾部)": ({"407": 'V'}, "Y407V", "Genentech", "【双抗组装】经典 KIH 凹陷端，扩大疏水口袋"),
    
    "EW 平台 (E)": ({"392": 'E'}, "K392E", "EW-RVT", "【双抗组装】电荷反转，与 Chain B 的 R 配对"),
    "EW 平台 (W)": ({"409": 'W'}, "K409W", "EW-RVT", "【双抗组装】引入疏水凸起，插入 Chain B 凹陷"),
    "RVT 平台 (R)": ({"347": 'R'}, "Q347R", "EW-RVT", "【双抗组装】引入正电荷，与 Chain A 的 E 配对"),
    "RVT 平台 (VT)": ({"399": 'V', "405": 'T'}, "D399V, F405T", "EW-RVT", "【双抗组装】构建疏水容纳凹陷"),

    "Azymetric Chain A (VY)": ({"350": 'V', "351": 'Y'}, "T350V, L351Y", "Zymeworks", "【双抗组装】不对称支架 A 链前段"),
    "Azymetric Chain A (AV)": ({"405": 'A', "407": 'V'}, "F405A, Y407V", "Zymeworks", "【双抗组装】不对称支架 A 链尾段"),
    "Azymetric Chain B (V)": ({"350": 'V'}, "T350V", "Zymeworks", "【双抗组装】不对称支架 B 链"),
    "Azymetric Chain B (L)": ({"366": 'L'}, "T366L", "Zymeworks", "【双抗组装】不对称支架 B 链核心"),
    "Azymetric Chain B (LW)": ({"392": 'L', "394": 'W'}, "K392L, T394W", "Zymeworks", "【双抗组装】不对称支架 B 链空间匹配"),

    "Charge Steer Chain A (K1)": ({"356": 'K'}, "E356K", "Chugai/XmAb", "【双抗组装】负转正电荷，利用静电排斥同源二聚"),
    "Charge Steer Chain A (K2)": ({"399": 'K'}, "D399K", "Chugai/XmAb", "【双抗组装】负转正电荷，深度静电互补"),
    "Charge Steer Chain B (D1)": ({"392": 'D'}, "K392D", "Chugai/XmAb", "【双抗组装】正转负电荷，配合 Chain A 的 K 突变"),
    "Charge Steer Chain B (D2)": ({"409": 'D'}, "K409D", "Chugai/XmAb", "【双抗组装】正转负电荷，配合 Chain A 的 K 突变"),

    "YTE": ({"252": 'Y', "254": 'T', "256": 'E'}, "M252Y, S254T, T256E", "AstraZeneca", "【PK优化】增强酸性下对 FcRn 的亲和力，延长半衰期"),
    "LS": ({"428": 'L', "434": 'S'}, "M428L, N434S", "Xencor", "【PK优化】当前最主流的长效修饰方案"),
    "N434A": ({"434": 'A'}, "N434A", "通用", "【PK优化】适度增强 FcRn 结合，延长半衰期"),
    "IHH (I253A)": ({"253": 'A'}, "I253A", "通用", "【快速清除】彻底破坏 FcRn 结合，极速缩短半衰期"),
    "IHH (H310A)": ({"310": 'A'}, "H310A", "通用", "【快速清除】协同阻止抗体被再循环"),

    "Protein A 破坏": ({"435": 'R', "436": 'F'}, "H435R, Y436F", "通用", "【CMC纯化】破坏与 Protein A 的结合，利用梯度洗脱纯化双抗"),
    "S228P": ({"228": 'P'}, "S228P", "通用", "【结构稳定】IgG4 专属，强制形成链间二硫键，防止 Fab 臂交换")
}

# ==========================================
# 3. Biopython 序列对齐与映射引擎
# ==========================================
def map_sequence_to_eu(target_seq):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.match_score = 2
    aligner.mismatch_score = -1

    alignments = aligner.align(WT_REF_SEQ, target_seq)
    if not alignments: return {}

    best_alignment = alignments[0]
    target_eu_map = {}
    
    aligned_ref = best_alignment.aligned[0]
    aligned_tgt = best_alignment.aligned[1]

    for (r_start, r_end), (t_start, t_end) in zip(aligned_ref, aligned_tgt):
        for i in range(r_end - r_start):
            eu_pos = str(WT_START_EU + r_start + i)
            target_eu_map[eu_pos] = target_seq[t_start + i]
            
    return target_eu_map

def check_mutation(pos_dict, eu_map):
    for pos, allowed_aa in pos_dict.items():
        if pos not in eu_map: return False
        if isinstance(allowed_aa, list):
            if eu_map[pos] not in allowed_aa: return False
        else:
            if eu_map[pos] != allowed_aa: return False
    return True

def analyze_fc(sequence):
    seq = sequence.upper().replace(" ", "").replace("\n", "")
    results = {"isotype": "未知 (未检测到标准骨架)", "mutations": [], "allotypes": [], "has_fc": False}
    
    eu_map = map_sequence_to_eu(seq)
    if not eu_map or len(eu_map) < 30: 
        return results
        
    results["has_fc"] = True

    for iso, motif in ISOTYPE_MOTIFS.items():
        if re.search(motif, seq):
            results["isotype"] = iso
            break
            
    for allo_name, (pos_req, eu_num, desc) in ALLOTYPE_DB.items():
        if check_mutation(pos_req, eu_map):
            results["allotypes"].append({"名称": allo_name, "位置": eu_num, "描述": desc})

    for mut_name, (pos_req, eu_num, platform, func) in MUTATION_DB.items():
        if check_mutation(pos_req, eu_map):
            results["mutations"].append({"突变简称": mut_name, "EU 编号": eu_num, "技术溯源": platform, "生物学功能": func})
            
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
raw_input = st.text_area("📥 粘贴抗体全长链或 Fc 段序列 (支持多条 FASTA，无惧移码/截断):", height=200)

if st.button("🔍 启动全境 Fc 深度解码", type="primary"):
    if raw_input:
        seq_dict = parse_fasta(raw_input)
        report_data = []
        deduction_reports = {} 
        
        with st.spinner("正在调用 Biopython 引擎比对绝对空间坐标并扫库..."):
            for name, seq in seq_dict.items():
                res = analyze_fc(seq)
                current_muts = [m["突变简称"] for m in res["mutations"]]
                current_allos = [a["名称"] for a in res["allotypes"]]
                
                deduction_reports[name] = {
                    "isotype": res["isotype"],
                    "muts_obj": res["mutations"], 
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
            st.session_state['fc_deduction'] = deduction_reports 
            # 重置之前的 3D 渲染状态，防止互相干扰
            st.session_state['render_3d_seq'] = None
            
            def highlight_rows(row):
                mut_str = str(row['特定突变识别'])
                if any(x in mut_str for x in ["Knob", "Hole", "EW", "RVT", "Azymetric", "Charge Steer"]): return ['background-color: #e3f2fd'] * len(row)
                if any(x in mut_str for x in ["LALA", "PAA", "P329G", "D265", "FEA", "Aglycosylation"]): return ['background-color: #fff3e0'] * len(row)
                if any(x in mut_str for x in ["GA-SD", "AL-IE", "HexaBody"]): return ['background-color: #ffebee'] * len(row)
                if any(x in mut_str for x in ["Protein A", "S228P"]): return ['background-color: #e8f5e9'] * len(row)
                if any(x in mut_str for x in ["YTE", "LS", "IHH", "N434A"]): return ['background-color: #f3e5f5'] * len(row)
                return [''] * len(row)

            st.markdown("### 📊 基于绝对映射的工程化分析报告")
            st.dataframe(df.style.apply(highlight_rows, axis=1), use_container_width=True)
            
            # --- 智能战略级推演 & 反向排雷 ---
            st.markdown("### 💡 独立分子战略推演与反向排雷预警")
            for seq_name, data in deduction_reports.items():
                if data["isotype"] == "未知 (未检测到标准骨架)": continue
                
                muts = data["muts"]
                allos = data["allos"]
                iso = data["isotype"]
                
                with st.expander(f"📌 情报解密: {seq_name}", expanded=True):
                    has_warning = False
                    
                    if iso == "IgG4 (野生型)" and "S228P" not in muts:
                        st.error("🚨 **反向排雷 [CMC风险]：缺失 S228P 稳定突变！** 检测到天然 IgG4 骨架，建议引入 S228P 避免 Fab 臂交换。")
                        has_warning = True
                        
                    is_bispecific = any(x in "".join(muts) for x in ["Knob", "Hole", "EW", "Azymetric", "Charge Steer"])
                    if is_bispecific and "Protein A 破坏" not in muts:
                        st.warning("⚠️ **反向排雷 [下游工艺风险]：缺失不对称纯化突变！** 建议增加 H435R/Y436F 以利用梯度洗脱纯化双抗。")
                        has_warning = True
                        
                    allo_str_joined = " ".join(allos)
                    if "G1m1" in allo_str_joined and "nG1m17" in allo_str_joined:
                        st.error("🚨 **反向排雷 [免疫原性风险]：同种异型冲突！** 非天然的人造嵌合体 (G1m1 + nG1m17)，极高 ADA 激发风险。")
                        has_warning = True
                        
                    if len(set(muts)) >= 4:
                        st.warning(f"⚠️ **反向排雷 [结构稳定性]：Fc 突变负荷过高 ({len(set(muts))}种)。** 极易导致 CH3 域微观折叠异常和 ADA 风险。")
                        has_warning = True

                    if not has_warning: st.success("✅ **排雷扫描通过**：结构稳健，未见明显的 CMC 缺陷或 Allotype 冲突。")

                    st.markdown("---")
                    st.markdown("##### 核心技术路线研判：")
                    if any("Azymetric" in m for m in muts): st.success("🎯 **双抗支架：Zymeworks (Azymetric) 空间支架。**")
                    elif any("Charge Steer" in m for m in muts): st.success("🎯 **双抗支架：静电反转驱动 (Chugai/XmAb等)。**")
                    elif any("EW" in m for m in muts) or any("RVT" in m for m in muts): st.success("🎯 **双抗支架：高阶 EW-RVT 平台。**")
                    elif any("Knob" in m for m in muts) or any("Hole" in m for m in muts): st.success("🎯 **双抗支架：经典 Knob-in-Hole (Genentech)。**")
                        
                    if any("LALA" in m for m in muts) and any("P329G" in m for m in muts): st.warning("🛡️ **杀伤机制：LALA-PG 终极效应子沉默。** (规避 CRS 的最佳选择)")
                    elif any("HexaBody" in m for m in muts): st.error("⚔️ **杀伤机制：补体风暴激发器 (Genmab HexaBody CDC)。**")
                    elif any("GA-SD" in m for m in muts) or any("AL-IE" in m for m in muts): st.error("⚔️ **杀伤机制：超级 ADCC 增强。**")
                        
                    if any("YTE" in m for m in muts) or any("LS" in m for m in muts): st.info("⏱️ **PK 设计：超长效修饰。** (大幅延长半衰期设计)")
                    elif any("IHH" in m for m in muts): st.info("☢️ **PK 设计：极速体内清除。** (通常用于核药/ADC)")
                    if not muts: st.markdown("- 🧬 **常规抗体**：未检测到特殊的工程化修饰意图。")
    else:
        st.error("请输入序列！")

# ==========================================
# 6. 🧊 3D 突变空间靶向映射实验室 (终极原生注入版)
# ==========================================
st.markdown("---")
st.markdown("### 🧊 3D 突变空间靶向映射实验室")

if 'fc_deduction' in st.session_state and st.session_state['fc_deduction']:
    st.info("💡 系统将加载标准人源 IgG1 Fc 的晶体结构 (PDB: 1FCC，A/B两链)，并将上方解析出的 EU 突变位点以【彩色球体】锚定在 3D 骨架上。")
    
    valid_seqs = [name for name, data in st.session_state['fc_deduction'].items() if data['muts_obj']]
    
    if valid_seqs:
        col_ctrl, col_graph = st.columns([1, 2])
        
        with col_ctrl:
            st.markdown("#### 控制面板")
            selected_seq = st.selectbox("🎯 选择候选链:", valid_seqs)
            
            # 按钮不再执行渲染，只负责将状态固化到 session_state 中
            if st.button("🏗️ 启动空间靶向渲染", type="primary", use_container_width=True):
                st.session_state['render_3d_seq'] = selected_seq
            
            st.markdown("""
            **颜色图例指引：**
            🔵 **蓝色球体**：异源二聚化骨架 (Knob, Hole 等)
            🟠 **橙色球体**：毒性沉默 (LALA, P329G, D265A 等)
            🔴 **红色球体**：杀伤功能增强 (HexaBody 等)
            🟢 **绿色球体**：CMC 稳定与纯化 (S228P, Protein A 破坏)
            🟣 **紫色球体**：PK 调节 (YTE, LS 等)
            """)
        
        with col_graph:
            # 只要 session_state 里有值，就跨越按钮生命周期，强制长久渲染！
            if st.session_state.get('render_3d_seq'):
                render_target = st.session_state['render_3d_seq']
                mut_data = st.session_state['fc_deduction'][render_target]['muts_obj']
                
                st.markdown(f"**🔬 当前空间靶向分子：** `{render_target}`")
                
                # 构建底层的 3D 引擎对象
                view = py3Dmol.view(width=800, height=500)
                view.addModel(query='pdb:1FCC')
                view.setStyle({'chain': 'A'}, {'cartoon': {'color': '#b0bec5'}})
                view.setStyle({'chain': 'B'}, {'cartoon': {'color': '#eceff1'}})
                view.setStyle({'chain': 'C'}, {'cartoon': {'hidden': True}})
                
                # 遍历打点
                for mut in mut_data:
                    mut_name = mut["突变简称"]
                    eu_str = mut["EU 编号"]
                    positions = re.findall(r'\d+', eu_str)
                    
                    color = 'red'
                    if any(x in mut_name for x in ["Knob", "Hole", "EW", "RVT", "Azymetric", "Charge Steer"]): color = '#2196f3'
                    elif any(x in mut_name for x in ["LALA", "PAA", "P329G", "D265", "FEA", "Aglycosylation"]): color = '#ff9800'
                    elif any(x in mut_name for x in ["GA-SD", "AL-IE", "HexaBody"]): color = '#e53935'
                    elif any(x in mut_name for x in ["Protein A", "S228P"]): color = '#4caf50'
                    elif any(x in mut_name for x in ["YTE", "LS", "IHH", "N434A"]): color = '#9c27b0'

                    for pos in positions:
                        view.addStyle({'resi': str(pos), 'chain': ['A', 'B']}, {'cartoon': {'color': '#cfd8dc'}, 'sphere': {'color': color, 'radius': 1.8}})
                        
                view.zoomTo()
                
                # 【终极奥义】：直接抽取 py3Dmol 的底层 HTML 代码，用 Streamlit 原生组件硬塞进网页！无视所有的 iframe 兼容 Bug！
                html_code = view._make_html()
                components.html(html_code, height=500, width=800)
    else:
        st.warning("⚠️ 刚才输入的序列中未检测到已知突变特征，无需空间映射。")
else:
    st.warning("⚠️ 请先在最上方输入序列文本并点击『启动全境 Fc 深度解码』。")
