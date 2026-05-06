"""
GeoCalc — Capacidad de Soporte + Asentamientos — Zapata Aislada
===============================================================
Capacidad de soporte: Terzaghi (1943), Meyerhof (1963), Hansen (1970), Vesić (1973)
Asentamientos: Elástico (inmediato), Consolidación primaria (Terzaghi),
               Consolidación secundaria, Schmertmann (arena)

Instalación:
    pip install streamlit plotly pandas

Correr:
    streamlit run app.py
"""

import math
import streamlit as st
import plotly.graph_objects as go
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURACIÓN DE PÁGINA
# ─────────────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="GeoCalc · Capacidad de Soporte + Asentamientos",
    page_icon="⬡",
    layout="wide",
)

st.markdown("""
<style>
    .titulo { font-size: 2rem; font-weight: 800; color: #e8c547; }
    .subtitulo { color: #6b7280; font-size: 0.85rem; margin-bottom: 24px; }
    .alerta-roja { background:#2d1515; border-left:4px solid #ff6b6b;
                   padding:12px 16px; border-radius:4px; color:#ff6b6b;
                   font-size:0.9rem; margin: 8px 0; }
    .alerta-verde { background:#152d1e; border-left:4px solid #4ecdc4;
                    padding:12px 16px; border-radius:4px; color:#4ecdc4;
                    font-size:0.9rem; margin: 8px 0; }
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="titulo">⬡ GeoCalc</div>', unsafe_allow_html=True)
st.markdown('<div class="subtitulo">Capacidad de Soporte + Asentamientos · Zapata Aislada</div>', unsafe_allow_html=True)
st.divider()


# ─────────────────────────────────────────────────────────────────────────────
# FUNCIONES — CAPACIDAD DE SOPORTE
# ─────────────────────────────────────────────────────────────────────────────

def factores_capacidad(phi_deg: float) -> dict:
    if phi_deg == 0:
        return {"Nc": 5.14, "Nq": 1.0, "Ng_M": 0.0, "Ng_H": 0.0, "Ng_V": 0.0}
    phi = math.radians(phi_deg)
    Nq   = math.exp(math.pi * math.tan(phi)) * math.tan(math.pi / 4 + phi / 2) ** 2
    Nc   = (Nq - 1) / math.tan(phi)
    Ng_M = (Nq - 1) * math.tan(math.radians(1.4 * phi_deg))
    Ng_H = 1.5 * (Nq - 1) * math.tan(phi)
    Ng_V = 2.0 * (Nq + 1) * math.tan(phi)
    return {"Nc": Nc, "Nq": Nq, "Ng_M": Ng_M, "Ng_H": Ng_H, "Ng_V": Ng_V}


def peso_efectivo(gamma, gamma2, nf, Df):
    if nf == "Sin nivel freático":
        return gamma, gamma
    elif nf == "En superficie":
        return gamma2, gamma2
    elif nf == "A nivel de desplante":
        return gamma, gamma2
    else:
        return gamma, (gamma + gamma2) / 2


def terzaghi_ult(tipo, B, L, Df, c, phi, gamma, gamma2, nf):
    f = factores_capacidad(phi)
    Nc, Nq, Ng = f["Nc"], f["Nq"], f["Ng_H"]
    q_g, b_g = peso_efectivo(gamma, gamma2, nf, Df)
    q = q_g * Df
    sc_map = {"cuadrada": 1.3, "circular": 1.3, "continua": 1.0, "rectangular": 1 + 0.2*(B/L)}
    sg_map = {"cuadrada": 0.8, "circular": 0.6, "continua": 1.0, "rectangular": max(1 - 0.2*(B/L), 0)}
    return max(0.0, c*Nc*sc_map[tipo] + q*Nq + 0.5*b_g*B*Ng*sg_map[tipo])


def meyerhof_ult(tipo, B, L, Df, c, phi, gamma, gamma2, nf):
    f = factores_capacidad(phi)
    Nc, Nq, Ng = f["Nc"], f["Nq"], f["Ng_M"]
    q_g, b_g = peso_efectivo(gamma, gamma2, nf, Df)
    q = q_g * Df
    phi_r = math.radians(phi)
    ratio = 1.0 if tipo == "circular" else (0.001 if tipo == "continua" else B/L)
    Fcs = 1 + ratio*(Nq/Nc) if Nc > 0 else 1
    Fqs = 1 + ratio*math.tan(phi_r)
    Fgs = max(1 - 0.4*ratio, 0)
    Dr = Df/B
    if phi == 0:
        Fcd = 1 + 0.4*Dr
        Fqd = 1.0
    else:
        Fqd = 1 + 2*math.tan(phi_r)*(1 - math.sin(phi_r))**2*Dr
        Fcd = Fqd - (1 - Fqd)/(Nc*math.tan(phi_r))
    return max(0.0, c*Nc*Fcs*Fcd + q*Nq*Fqs*Fqd + 0.5*b_g*B*Ng*Fgs)


def hansen_ult(tipo, B, L, Df, c, phi, gamma, gamma2, nf):
    f = factores_capacidad(phi)
    Nc, Nq, Ng = f["Nc"], f["Nq"], f["Ng_H"]
    q_g, b_g = peso_efectivo(gamma, gamma2, nf, Df)
    q = q_g * Df
    phi_r = math.radians(phi)
    ratio = 1.0 if tipo == "circular" else (0.001 if tipo == "continua" else B/L)
    sq = 1 + ratio*math.sin(phi_r)
    sc = (1 + 0.4*ratio) if phi == 0 else (1 + ratio*(Nq/Nc) if Nc > 0 else 1)
    sg = max(1 - 0.4*ratio, 0)
    Dr = Df/B
    atan_Dr = math.atan(Dr)
    if phi == 0:
        dc = 1 + 0.4*atan_Dr
        dq = 1.0
    else:
        dq = 1 + 2*math.tan(phi_r)*(1 - math.sin(phi_r))**2*atan_Dr
        dc = dq - (1 - dq)/(Nc*math.tan(phi_r))
    return max(0.0, c*Nc*sc*dc + q*Nq*sq*dq + 0.5*b_g*B*Ng*sg)


def vesic_ult(tipo, B, L, Df, c, phi, gamma, gamma2, nf):
    f = factores_capacidad(phi)
    Nc, Nq, Ng = f["Nc"], f["Nq"], f["Ng_V"]
    q_g, b_g = peso_efectivo(gamma, gamma2, nf, Df)
    q = q_g * Df
    phi_r = math.radians(phi)
    ratio = 1.0 if tipo == "circular" else (0.001 if tipo == "continua" else B/L)
    Fcs = 1 + ratio*(Nq/Nc) if Nc > 0 else 1
    Fqs = 1 + ratio*math.tan(phi_r)
    Fgs = max(1 - 0.4*ratio, 0)
    Dr = Df/B
    atan_Dr = math.atan(Dr)
    if phi == 0:
        Fcd = 1 + 0.4*atan_Dr
        Fqd = 1.0
    else:
        Fqd = 1 + 2*math.tan(phi_r)*(1 - math.sin(phi_r))**2*atan_Dr
        Fcd = Fqd - (1 - Fqd)/(Nc*math.tan(phi_r))
    return max(0.0, c*Nc*Fcs*Fcd + q*Nq*Fqs*Fqd + 0.5*b_g*B*Ng*Fgs)


# ─────────────────────────────────────────────────────────────────────────────
# FUNCIONES — ASENTAMIENTOS
# ─────────────────────────────────────────────────────────────────────────────

def asentamiento_elastico(q_neta, B, L, Es, nu, Df, tipo_zapata) -> float:
    """
    Asentamiento elástico (inmediato) — Método de Bowles / Das
    Se = q_neta * B * (1 - nu^2) / Es * F1 * Fox
    Resultado en mm.
    """
    if Es <= 0:
        return 0.0
    factores_forma = {
        "cuadrada":    0.82,
        "circular":    0.85,
        "continua":    1.52,
        "rectangular": 1.06,
    }
    F1  = factores_forma.get(tipo_zapata, 1.0)
    Fox = max(1 - 0.4 * math.log10(1 + Df / B), 0.5) if B > 0 else 1.0
    Se  = q_neta * B * (1 - nu**2) / Es * F1 * Fox
    return Se * 1000  # mm


def asentamiento_consolidacion_primaria(q_neta, H, Cc, e0, sigma_p, sigma_v0) -> float:
    """
    Consolidación primaria — Terzaghi
    Maneja suelos NC y OC.
    Resultado en mm.
    """
    if Cc <= 0 or H <= 0 or sigma_v0 <= 0:
        return 0.0
    delta_sigma = q_neta
    sigma_final = sigma_v0 + delta_sigma
    if sigma_v0 >= sigma_p:
        # Normalmente consolidado
        Sc = (Cc / (1 + e0)) * H * math.log10(sigma_final / sigma_v0)
    else:
        Cs = Cc / 5
        if sigma_final <= sigma_p:
            Sc = (Cs / (1 + e0)) * H * math.log10(sigma_final / sigma_v0)
        else:
            Sc = (Cs / (1 + e0)) * H * math.log10(sigma_p / sigma_v0) + \
                 (Cc / (1 + e0)) * H * math.log10(sigma_final / sigma_p)
    return max(Sc, 0.0) * 1000  # mm


def asentamiento_consolidacion_secundaria(Calpha, e0, H, t1, t2) -> float:
    """
    Consolidación secundaria (creep)
    Ss = (Calpha / (1+e0)) * H * log10(t2/t1)
    Resultado en mm.
    """
    if t2 <= t1 or H <= 0:
        return 0.0
    Calpha_prima = Calpha / (1 + e0)
    Ss = Calpha_prima * H * math.log10(t2 / t1)
    return max(Ss, 0.0) * 1000  # mm


def asentamiento_schmertmann(q_neta, B, Es, Df, gamma, tipo_zapata) -> float:
    """
    Método de Schmertmann (1970) — para arenas
    Se = C1 * C2 * q_neta * sum(Iz/Es * dz)
    Resultado en mm.
    """
    if Es <= 0 or q_neta <= 0:
        return 0.0
    sigma_v0 = gamma * Df
    C1 = max(1 - 0.5 * (sigma_v0 / q_neta), 0.5)
    C2 = 1 + 0.2 * math.log10(10 / 0.1)  # t = 10 años
    z_max = 2 * B if tipo_zapata in ["cuadrada", "circular"] else 4 * B
    n_capas = 20
    dz = z_max / n_capas
    suma = 0.0
    for i in range(n_capas):
        z_mid  = (i + 0.5) * dz
        z_peak = z_max / 4
        if z_mid <= z_peak:
            Iz = 0.6 * (z_mid / z_peak)
        else:
            Iz = 0.6 * (1 - (z_mid - z_peak) / (z_max - z_peak))
        Iz = max(Iz, 0)
        suma += (Iz / Es) * dz
    Se = C1 * C2 * q_neta * suma
    return max(Se, 0.0) * 1000  # mm


# ─────────────────────────────────────────────────────────────────────────────
# PANEL LATERAL — PARÁMETROS
# ─────────────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("📐 Geometría")

    tipo_zapata = st.selectbox(
        "Tipo de zapata",
        ["cuadrada", "rectangular", "circular", "continua"],
        format_func=lambda x: x.capitalize()
    )
    B  = st.number_input("B — Ancho (m)", value=1.5, min_value=0.1, step=0.1, format="%.2f")
    L  = st.number_input("L — Largo (m)", value=2.0, min_value=0.1, step=0.1, format="%.2f") \
         if tipo_zapata == "rectangular" else B
    Df = st.number_input("Df — Prof. desplante (m)", value=1.2, min_value=0.1, step=0.1, format="%.2f")

    st.divider()
    st.header("🪨 Suelo — Resistencia")

    tipo_suelo = st.selectbox(
        "Tipo de suelo",
        ["Mixto (c y φ)", "Cohesivo puro (arcilla φ=0)", "Friccionante (arena c=0)"]
    )
    if tipo_suelo == "Cohesivo puro (arcilla φ=0)":
        c   = st.number_input("c — Cohesión (kPa)", value=50.0, min_value=0.0, step=1.0)
        phi = 0
        st.info("φ = 0°")
    elif tipo_suelo == "Friccionante (arena c=0)":
        c   = 0.0
        phi = st.slider("φ — Fricción (°)", 15, 45, 30)
        st.info("c = 0 kPa")
    else:
        c   = st.number_input("c — Cohesión (kPa)", value=20.0, min_value=0.0, step=1.0)
        phi = st.slider("φ — Fricción (°)", 0, 45, 25)

    gamma = st.number_input("γ — Peso unitario (kN/m³)", value=18.0, min_value=10.0, max_value=25.0, step=0.5)

    nf = st.selectbox("Nivel freático",
        ["Sin nivel freático", "En superficie", "A nivel de desplante", "A B metros bajo desplante"])
    gamma2 = st.number_input("γ' — Sumergido (kN/m³)", value=9.0, min_value=6.0, max_value=14.0, step=0.5) \
             if nf != "Sin nivel freático" else gamma

    FS = st.slider("Factor de seguridad (FS)", 1.5, 4.0, 3.0, 0.5)

    st.divider()
    st.header("⚡ Carga aplicada")
    q_aplicada = st.number_input(
        "q — Carga neta de la zapata (kPa)",
        value=150.0, min_value=1.0, step=10.0,
        help="Carga neta = Carga total / Área − γ·Df"
    )

    st.divider()
    st.header("📉 Suelo — Asentamientos")

    Es = st.number_input(
        "Es — Módulo de elasticidad (kPa)",
        value=15000.0, min_value=500.0, step=500.0,
        help="Arena suelta: 5,000–15,000 | Arena densa: 20,000–80,000 | Arcilla: 2,000–20,000"
    )
    nu = st.slider("ν — Coeficiente de Poisson", 0.1, 0.5, 0.3, 0.05,
                   help="Arena: 0.25–0.35 | Arcilla saturada: 0.45–0.5")

    st.markdown("**Consolidación primaria**")
    Cc = st.number_input("Cc — Índice de compresión", value=0.35, min_value=0.01, step=0.01, format="%.3f")
    e0 = st.number_input("e0 — Relación de vacíos", value=0.80, min_value=0.1, step=0.05, format="%.3f")
    H  = st.number_input("H — Espesor capa compresible (m)", value=4.0, min_value=0.1, step=0.5)
    sigma_v0 = st.number_input("σ'v0 — Esfuerzo efectivo inicial (kPa)", value=60.0, min_value=1.0, step=5.0,
                                help="γ' × z al centro de la capa compresible")
    sigma_p  = st.number_input("σ'p — Presión de preconsolidación (kPa)", value=80.0, min_value=1.0, step=5.0,
                                help="Si σ'p ≤ σ'v0 → NC | Si σ'p > σ'v0 → OC")

    st.markdown("**Consolidación secundaria**")
    Calpha = st.number_input("Cα — Coef. secundario", value=0.02, min_value=0.001, step=0.005, format="%.4f",
                              help="Turba: 0.05–0.1 | Arcilla blanda: 0.01–0.03 | Arcilla firme: 0.005–0.01")
    t1 = st.number_input("t₁ — Fin consol. primaria (años)", value=1.0, min_value=0.1, step=0.5)
    t2 = st.number_input("t₂ — Tiempo de análisis (años)", value=30.0, min_value=1.0, step=1.0)


# ─────────────────────────────────────────────────────────────────────────────
# CÁLCULOS
# ─────────────────────────────────────────────────────────────────────────────
args_cap = (tipo_zapata, B, L, Df, c, phi, gamma, gamma2, nf)

qu = {
    "Terzaghi": terzaghi_ult(*args_cap),
    "Meyerhof": meyerhof_ult(*args_cap),
    "Hansen":   hansen_ult(*args_cap),
    "Vesić":    vesic_ult(*args_cap),
}
qa = {k: v/FS for k, v in qu.items()}

colores = {
    "Terzaghi": "#e8c547",
    "Meyerhof": "#4ecdc4",
    "Hansen":   "#ff6b6b",
    "Vesić":    "#a78bfa",
}

valores_qu  = list(qu.values())
promedio_qu = sum(valores_qu) / len(valores_qu)
promedio_qa = promedio_qu / FS
min_qu = min(valores_qu)
max_qu = max(valores_qu)
f = factores_capacidad(phi)

# Asentamientos
Se  = asentamiento_elastico(q_aplicada, B, L, Es, nu, Df, tipo_zapata)
Sc  = asentamiento_consolidacion_primaria(q_aplicada, H, Cc, e0, sigma_p, sigma_v0)
Ss  = asentamiento_consolidacion_secundaria(Calpha, e0, H, t1, t2)
Sch = asentamiento_schmertmann(q_aplicada, B, Es, Df, gamma, tipo_zapata)
St  = Se + Sc + Ss

LIMITE_TOTAL = 25.0  # mm


# ─────────────────────────────────────────────────────────────────────────────
# TABS PRINCIPALES
# ─────────────────────────────────────────────────────────────────────────────
tab1, tab2 = st.tabs(["⬡ Capacidad de Soporte", "📉 Asentamientos"])


# ═══════════════════════════════════════════════════════════════════════════
# TAB 1 — CAPACIDAD DE SOPORTE
# ═══════════════════════════════════════════════════════════════════════════
with tab1:

    st.subheader("Resumen")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("qu promedio", f"{promedio_qu:.1f} kPa")
    c2.metric("qa promedio", f"{promedio_qa:.1f} kPa")
    c3.metric("qu mínimo",   f"{min_qu:.1f} kPa")
    c4.metric("qu máximo",   f"{max_qu:.1f} kPa")

    if q_aplicada <= promedio_qa:
        st.markdown(f'<div class="alerta-verde">✅ Carga aplicada ({q_aplicada:.1f} kPa) ≤ qa promedio ({promedio_qa:.1f} kPa) — Diseño SEGURO</div>', unsafe_allow_html=True)
    else:
        st.markdown(f'<div class="alerta-roja">⚠️ Carga aplicada ({q_aplicada:.1f} kPa) > qa promedio ({promedio_qa:.1f} kPa) — Revisar diseño</div>', unsafe_allow_html=True)

    st.divider()
    st.subheader("Resultados por Método")

    notas = {
        "Terzaghi": "Clásico y conservador",
        "Meyerhof": "Forma + profundidad",
        "Hansen":   "Recomendado arcillas",
        "Vesić":    "Nγ mayor → conservador",
    }
    cols = st.columns(4)
    for col, (nombre, qu_val) in zip(cols, qu.items()):
        with col:
            tag = "🔻 MÍN" if qu_val == min_qu else ("🔺 MÁX" if qu_val == max_qu else "")
            st.markdown(f"**{nombre}** {tag}")
            st.markdown(f"<span style='color:{colores[nombre]};font-size:1.6rem;font-weight:800'>{qu_val:.1f}</span> kPa", unsafe_allow_html=True)
            st.caption("qu — capacidad última")
            st.metric(label=f"qa (÷FS {FS})", value=f"{qa[nombre]:.1f} kPa")
            st.caption(notas[nombre])

    st.divider()
    st.subheader("Comparación Visual")

    fig1 = go.Figure(go.Bar(
        x=list(qu.keys()),
        y=list(qu.values()),
        marker_color=list(colores.values()),
        text=[f"{v:.1f}" for v in qu.values()],
        textposition="outside",
    ))
    fig1.add_hline(y=promedio_qu, line_dash="dash",
                   line_color="rgba(255,255,255,0.3)",
                   annotation_text=f"Promedio: {promedio_qu:.1f} kPa",
                   annotation_font_color="#9ca3af")
    fig1.add_hline(y=q_aplicada, line_dash="dot",
                   line_color="rgba(255,107,107,0.6)",
                   annotation_text=f"Carga aplicada: {q_aplicada:.1f} kPa",
                   annotation_font_color="#ff6b6b")
    fig1.update_layout(
        plot_bgcolor="#13161d", paper_bgcolor="#13161d",
        font_color="#e8eaf0", yaxis_title="qu (kPa)",
        showlegend=False, height=380,
        margin=dict(t=40, b=20),
        yaxis=dict(gridcolor="#252a38"),
    )
    st.plotly_chart(fig1, use_container_width=True)

    st.subheader("Tabla de Resultados")
    df_cap = pd.DataFrame({
        "Método":   list(qu.keys()),
        "qu (kPa)": [round(v, 2) for v in qu.values()],
        "qa (kPa)": [round(v, 2) for v in qa.values()],
        "FS real":  [round(v/q_aplicada, 2) if q_aplicada > 0 else 0 for v in qu.values()],
        "Nc": [round(f["Nc"], 3)]*4,
        "Nq": [round(f["Nq"], 3)]*4,
        "Nγ": [round(f["Ng_M"],3), round(f["Ng_M"],3), round(f["Ng_H"],3), round(f["Ng_V"],3)],
    })
    st.dataframe(df_cap, hide_index=True, use_container_width=True)
    st.download_button("⬇ Descargar CSV — Capacidad",
                       df_cap.to_csv(index=False).encode("utf-8"),
                       "capacidad_soporte.csv", "text/csv")


# ═══════════════════════════════════════════════════════════════════════════
# TAB 2 — ASENTAMIENTOS
# ═══════════════════════════════════════════════════════════════════════════
with tab2:

    st.subheader("Resumen de Asentamientos")

    ca1, ca2, ca3, ca4, ca5 = st.columns(5)
    ca1.metric("Se — Elástico",           f"{Se:.1f} mm")
    ca2.metric("Sc — Consol. primaria",   f"{Sc:.1f} mm")
    ca3.metric("Ss — Consol. secundaria", f"{Ss:.1f} mm")
    ca4.metric("Sch — Schmertmann",       f"{Sch:.1f} mm")
    ca5.metric("St — Total (Se+Sc+Ss)",   f"{St:.1f} mm")

    if St <= LIMITE_TOTAL:
        st.markdown(f'<div class="alerta-verde">✅ Asentamiento total ({St:.1f} mm) dentro del límite típico de {LIMITE_TOTAL} mm</div>', unsafe_allow_html=True)
    else:
        st.markdown(f'<div class="alerta-roja">⚠️ Asentamiento total ({St:.1f} mm) supera el límite de {LIMITE_TOTAL} mm — Revisar diseño</div>', unsafe_allow_html=True)

    st.caption("Referencia: asentamiento total ≤ 25 mm | diferencial ≤ 19 mm (edificaciones típicas, Das 2011)")

    st.divider()
    st.subheader("Comparación por Componente")

    fig2 = go.Figure(go.Bar(
        x=["Elástico (Se)", "Consol. Primaria (Sc)", "Consol. Secundaria (Ss)", "Schmertmann (Sch)"],
        y=[Se, Sc, Ss, Sch],
        marker_color=["#e8c547", "#4ecdc4", "#ff6b6b", "#a78bfa"],
        text=[f"{v:.1f} mm" for v in [Se, Sc, Ss, Sch]],
        textposition="outside",
    ))
    fig2.add_hline(y=LIMITE_TOTAL, line_dash="dash",
                   line_color="rgba(255,107,107,0.5)",
                   annotation_text=f"Límite típico: {LIMITE_TOTAL} mm",
                   annotation_font_color="#ff6b6b")
    fig2.update_layout(
        plot_bgcolor="#13161d", paper_bgcolor="#13161d",
        font_color="#e8eaf0", yaxis_title="Asentamiento (mm)",
        showlegend=False, height=380,
        margin=dict(t=40, b=20),
        yaxis=dict(gridcolor="#252a38"),
    )
    st.plotly_chart(fig2, use_container_width=True)

    st.divider()
    st.subheader("Evolución del Asentamiento en el Tiempo")

    tiempos = [0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 30, 50, 100]
    asen_tiempo = []
    for t in tiempos:
        if t <= t1:
            fraccion = t / t1
            asen_tiempo.append(Se + Sc * fraccion)
        else:
            Ss_t = asentamiento_consolidacion_secundaria(Calpha, e0, H, t1, t)
            asen_tiempo.append(Se + Sc + Ss_t)

    fig3 = go.Figure()
    fig3.add_trace(go.Scatter(
        x=tiempos, y=asen_tiempo,
        mode="lines+markers",
        line=dict(color="#e8c547", width=2),
        marker=dict(size=6),
        name="Asentamiento total"
    ))
    fig3.add_hline(y=LIMITE_TOTAL, line_dash="dash",
                   line_color="rgba(255,107,107,0.5)",
                   annotation_text="Límite 25 mm",
                   annotation_font_color="#ff6b6b")
    fig3.add_vline(x=t1, line_dash="dot",
                   line_color="rgba(78,205,196,0.5)",
                   annotation_text=f"Fin consol. primaria (t₁={t1} a)",
                   annotation_font_color="#4ecdc4")
    fig3.update_layout(
        plot_bgcolor="#13161d", paper_bgcolor="#13161d",
        font_color="#e8eaf0",
        xaxis_title="Tiempo (años)", yaxis_title="Asentamiento (mm)",
        height=400, margin=dict(t=20, b=20),
        xaxis=dict(gridcolor="#252a38", type="log"),
        yaxis=dict(gridcolor="#252a38"),
    )
    st.plotly_chart(fig3, use_container_width=True)

    st.divider()
    st.subheader("Tabla de Resultados — Asentamientos")

    estado_suelo = "NC (normalmente consolidado)" if sigma_v0 >= sigma_p else "OC (preconsolidado)"
    df_asen = pd.DataFrame({
        "Componente":  ["Elástico (Se)", "Consol. Primaria (Sc)", "Consol. Secundaria (Ss)", "Schmertmann (Sch)", "Total Se+Sc+Ss"],
        "Valor (mm)":  [round(Se,2), round(Sc,2), round(Ss,2), round(Sch,2), round(St,2)],
        "Método":      ["Bowles/Das", "Terzaghi", "Creep", "Schmertmann (1970)", "—"],
        "Aplicación":  ["General", estado_suelo, "General", "Arenas", "—"],
    })
    st.dataframe(df_asen, hide_index=True, use_container_width=True)

    with st.expander("📘 Estado de consolidación del suelo"):
        st.markdown(f"""
        - **σ'v0** = {sigma_v0:.1f} kPa &nbsp;|&nbsp; **σ'p** = {sigma_p:.1f} kPa
        - **Estado**: {estado_suelo}
        - Cs estimado (índice de recompresión) = **{round(Cc/5,4)}** (= Cc/5)
        - Si **σ'v0 ≥ σ'p** → suelo NC → se usa Cc completo
        - Si **σ'v0 < σ'p** → suelo OC → se usa Cs en zona elástica y Cc al superar σ'p
        """)

    st.download_button("⬇ Descargar CSV — Asentamientos",
                       df_asen.to_csv(index=False).encode("utf-8"),
                       "asentamientos.csv", "text/csv")

st.divider()
st.caption("GeoCalc · Das, B.M. — Principles of Foundation Engineering · Bowles, J.E. — Foundation Analysis and Design · Uso académico y profesional")
