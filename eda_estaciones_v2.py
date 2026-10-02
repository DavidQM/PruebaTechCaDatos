#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EDA y detección de anomalías de las estaciones 35, 367, 803 y 6004
===================================================================

Tarea A  Análisis exploratorio (A1 a A5): tablas y gráficas por estación y comparativas.
Tarea C  Detección de anomalías:
    C1  reglas deterministas (rango, salto, persistencia, coherencia)
    C2  métodos estadísticos y de ML (residuos con MAD, Isolation Forest, cambios de régimen)
    C3  comparación entre enfoques y contra la columna `calidad` (etiqueta débil)
    C4  clasificación de cada anomalía relevante
    C5  esquema de banderas y diccionario de los códigos de calidad inferidos

Uso
    python eda_anomalias_estaciones.py --datos RUTA_CSV --salida resultados
    python eda_anomalias_estaciones.py --datos RUTA_CSV --solo-eda
    python eda_anomalias_estaciones.py --datos RUTA_CSV --solo-anomalias --lat 6.25 --lon -75.57

Dependencias
    pip install numpy pandas scipy scikit-learn matplotlib

Entradas
    Cuatro CSV con columnas codigo, fecha_hora, <variables>, calidad, a 1 minuto:
        *pluviometrica_35_*.csv   (p1, p2)
        *meteorologica_367_*.csv  (h, t, pr, vv, vv_max, dv, dv_max, p)
        *nivel_803_*.csv          (nivel)
        *piranometro_6004_*.csv   (radiacion)

Salidas (carpeta --salida)
    fig_*.png     gráficas de la Tarea A y de la Tarea C
    tabla_*.csv   tablas de resultados
    datos_banderas_<estacion>.csv.gz   datos con las banderas asignadas (con --exportar-datos)

Todo umbral que no proviene de una norma está en los diccionarios de configuración
de la sección 1, para ajustarlo a la estación.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import warnings
from collections import OrderedDict

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from scipy import ndimage, stats
from scipy.signal import find_peaks
from sklearn.ensemble import IsolationForest
from sklearn.metrics import cohen_kappa_score
from sklearn.preprocessing import RobustScaler

warnings.filterwarnings("ignore")

# =====================================================================================
# 1. CONFIGURACIÓN
# =====================================================================================
ESTACIONES = OrderedDict([
    ("35", dict(nombre="Pluviométrica 35", vars=["p1", "p2"], principal="p1")),
    ("367", dict(nombre="Meteorológica 367", vars=["h", "t", "pr", "vv", "vv_max", "dv", "dv_max", "p"], principal="t")),
    ("803", dict(nombre="Nivel 803", vars=["nivel"], principal="nivel")),
    ("6004", dict(nombre="Piranómetro 6004", vars=["radiacion"], principal="radiacion")),
])
COLOR = {"35": "#2a78d6", "367": "#eb6834", "803": "#1baf7a", "6004": "#6250d6"}
GRIS, ROJO = "#6b7280", "#e34948"
CENTINELA = -999

# Rango físico (mínimo, máximo). None = sin cota. El nivel no tiene unidad declarada.
LIMITES = {
    "35": {"p1": (0, 5), "p2": (0, 5)},                                  # mm por minuto
    "367": {"h": (0, 100), "t": (-5, 45), "pr": (750, 900), "vv": (0, 60), "vv_max": (0, 80),
            "dv": (0, 360), "dv_max": (0, 360), "p": (0, 5)},            # pr: ajustar a la altitud
    "803": {"nivel": (0, 500)},
    "6004": {"radiacion": (0, None)},                                    # el máximo sale de la geometría solar
}
# Salto máximo entre minutos consecutivos (unidades de la variable por minuto)
SALTO_MAX = {"367": {"t": 2.0, "h": 20.0, "pr": 1.5, "vv": 10.0}, "803": {"nivel": 20.0}, "6004": {"radiacion": 800.0}}
SALTO_HORARIO_T = 3.0            # °C entre medias horarias consecutivas (criterio de los manuales del IDEAM)
# Persistencia: minutos consecutivos con el mismo valor
CONGELADO_MIN = {"367": {"t": 60, "h": 60, "pr": 60, "vv": 60, "dv": 60}, "803": {"nivel": 30}, "6004": {"radiacion": 10}}

P = dict(
    lat=6.25, lon=-75.57, meridiano=-75.0,   # supuestos de ubicación (Valle de Aburrá), parametrizables
    kt_max=1.20, noche_wm2=5.0,              # radiación: cota sobre la extraterrestre y umbral nocturno
    humedad100_min=30, lluvia_minima_min=15, # persistencia específica
    discrep_pulsos=2, cociente_tol=0.10, cociente_min_mm=2.0,
    nivel_alto_z=10.0, lluvia_previa_h=6, lluvia_previa_mm=1.0,
    mad_k=6.0, mad_ventana=181,
    if_contaminacion=0.005, if_arboles=200, if_max_muestras=20000, if_semilla=42,
    regimen_ventana_h=12, regimen_k=6.0,
    brecha_eventos=10, tolerancia_min=2, codigos_ok=(1,),
)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.titlesize": 9, "axes.titleweight": "bold",
                     "axes.titlelocation": "left", "axes.grid": True, "grid.alpha": .25,
                     "grid.linewidth": .6, "figure.dpi": 160})


# =====================================================================================
# 2. CARGA DE DATOS Y UTILIDADES
# =====================================================================================
def irradiancia_extraterrestre(idx_local, lat, lon, meridiano=-75.0):
    """Irradiancia horizontal fuera de la atmósfera (W/m²) para un índice en hora local."""
    doy = idx_local.dayofyear.values
    hora = idx_local.hour.values + idx_local.minute.values / 60
    b = 2 * np.pi * (doy - 81) / 364
    eot = 9.87 * np.sin(2 * b) - 7.53 * np.cos(b) - 1.5 * np.sin(b)           # ecuación del tiempo (min)
    hora_solar = hora + (lon - meridiano) / 15 + eot / 60
    decl = np.deg2rad(23.45) * np.sin(2 * np.pi * (284 + doy) / 365)
    ha = np.deg2rad(15 * (hora_solar - 12))
    la = np.deg2rad(lat)
    cosz = np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha)
    g0 = 1361 * (1 + 0.033 * np.cos(2 * np.pi * doy / 365)) * np.clip(cosz, 0, None)
    return pd.Series(g0, index=idx_local)


def estimar_desfase_horario(serie, lat, lon, rango=range(-12, 13)):
    """Desfase (h) que maximiza la correlación entre la radiación y la geometría solar."""
    base = serie.dropna()
    corr = {}
    for s in rango:
        g = irradiancia_extraterrestre(base.index + pd.Timedelta(hours=s), lat, lon)
        corr[s] = float(np.corrcoef(base.values, g.values)[0, 1])
    return max(corr, key=corr.get), corr


def racha_mask(s, n):
    """True en los minutos que pertenecen a una racha de >= n valores idénticos (NaN corta la racha)."""
    s = s.round(4)
    cambio = (s != s.shift()) | s.isna()
    gid = cambio.cumsum()
    tam = s.groupby(gid).transform("size")
    return (tam >= n) & s.notna()


def tramos(mascara, brecha=0):
    """Pares (inicio, fin) de posiciones de runs True; une los separados por <= brecha minutos."""
    m = np.asarray(mascara, bool)
    pos = np.flatnonzero(m)
    if pos.size == 0:
        return np.empty((0, 2), int)
    cortes = np.flatnonzero(np.diff(pos) > brecha + 1)
    return np.c_[np.r_[pos[0], pos[cortes + 1]], np.r_[pos[cortes], pos[-1]]]


def dilatar(mascara, k):
    return ndimage.binary_dilation(np.asarray(mascara, bool), structure=np.ones(2 * k + 1, bool))


class Datos:
    """Carga los CSV y los deja reindexados a una malla regular de 1 minuto."""

    def __init__(self, carpeta, desfase_6004="auto"):
        self.carpeta = carpeta
        crudo = {}
        for k in ESTACIONES:
            rutas = sorted(glob.glob(os.path.join(carpeta, f"*_{k}_*.csv")))
            if not rutas:
                raise FileNotFoundError(f"No se encontró el CSV de la estación {k} en {carpeta}")
            r = pd.read_csv(rutas[0])
            r["fecha_hora"] = pd.to_datetime(r["fecha_hora"])
            crudo[k] = r
        self.crudo = crudo
        ini = min(r.fecha_hora.min() for r in crudo.values()).floor("D")
        fin = max(r.fecha_hora.max() for r in crudo.values()).floor("D") + pd.Timedelta(days=1)
        self.idx = pd.date_range(ini, fin - pd.Timedelta(minutes=1), freq="1min")
        self.ini, self.fin = ini, fin
        self.df, self.sent, self.dup = {}, {}, {}
        for k, cfg in ESTACIONES.items():
            r = crudo[k].copy()
            r["_min"] = r.fecha_hora.dt.floor("1min")
            r["_irreg"] = r.fecha_hora.dt.second != 0
            dup = r["_min"].duplicated(keep=False)
            self.dup[k] = pd.Series(r.loc[dup, "_min"].unique(), dtype="datetime64[ns]")
            r = r.drop_duplicates("_min", keep="first").set_index("_min")
            d = r[cfg["vars"] + ["calidad"]].reindex(self.idx)
            d["presente"] = self.idx.isin(r.index)
            d["seg_irregular"] = r["_irreg"].reindex(self.idx).fillna(False).astype(bool)
            self.df[k] = d
            self.sent[k] = (d[cfg["vars"]] == CENTINELA).any(axis=1)
            self.dup[k] = pd.Series(self.idx.isin(self.dup[k]), index=self.idx)
        # desfase horario del piranómetro y cota solar
        rad = self.vals("6004")["radiacion"]
        if str(desfase_6004).lower() == "auto":
            self.desfase_6004, self.corr_desfase = estimar_desfase_horario(rad, P["lat"], P["lon"])
        else:
            self.desfase_6004 = int(desfase_6004)
            _, self.corr_desfase = estimar_desfase_horario(rad, P["lat"], P["lon"])
        g = irradiancia_extraterrestre(self.idx + pd.Timedelta(hours=self.desfase_6004), P["lat"], P["lon"])
        g.index = self.idx
        self.g0 = g

    def vals(self, k):
        """Variables con el centinela convertido a NaN, sobre la malla completa."""
        d = self.df[k][ESTACIONES[k]["vars"]]
        return d.where(d != CENTINELA)

    def etiqueta_debil(self, k):
        """True cuando la columna `calidad` trae un código distinto de los aceptados (solo filas presentes)."""
        c = self.df[k]["calidad"]
        return c.notna() & ~c.isin(P["codigos_ok"])

    def faltantes(self, k):
        """Minuto sin dato útil: sin registro o con el centinela."""
        return ~self.df[k]["presente"] | self.sent[k]


def serie_diaria(D, k):
    d = D.vals(k)
    if k == "35":
        return d.p1.resample("D").sum(min_count=1)
    if k == "803":
        return d.nivel.resample("D").median()
    if k == "6004":
        x = d.radiacion.resample("D").mean()
        return x.where(d.radiacion.notna().resample("D").mean() >= 0.9)
    return d.t.resample("D").mean()


def ciclo_diurno(D, k, desfase_h=0):
    s = D.vals(k)[ESTACIONES[k]["principal"]].copy()
    s.index = s.index + pd.Timedelta(hours=desfase_h)
    g = s.groupby(s.index.hour)
    return pd.DataFrame({"p25": g.quantile(.25), "p50": g.median(), "p75": g.quantile(.75), "media": g.mean()})


# =====================================================================================
# 3. TAREA A: EDA
# =====================================================================================
def tablas_eda(D, out):
    perfil, huecos, codigos, tend = [], [], [], []
    for k, cfg in ESTACIONES.items():
        raw, d = D.crudo[k], D.df[k]
        dif = raw.fecha_hora.diff().dropna().dt.total_seconds()
        for v in cfg["vars"]:
            x = raw[v]
            xv = x[x != CENTINELA]
            perfil.append(dict(estacion=k, variable=v, tipo=str(x.dtype), filas=len(raw), esperadas=len(D.idx),
                               orden_creciente=bool(raw.fecha_hora.is_monotonic_increasing),
                               ts_duplicados=int(raw.fecha_hora.duplicated().sum()),
                               intervalos_distintos_60s=int((dif != 60).sum()),
                               segundos_no_cero=int((raw.fecha_hora.dt.second != 0).sum()),
                               nulos=int(x.isna().sum()), centinela=int((x == CENTINELA).sum()),
                               negativos=int((xv < 0).sum()), ceros=int((xv == 0).sum()),
                               minimo=xv.min(), p01=xv.quantile(.01), mediana=xv.median(), p99=xv.quantile(.99),
                               maximo=xv.max(), media=xv.mean(), desviacion=xv.std(), asimetria=stats.skew(xv.dropna())))
        miss = D.faltantes(k)
        g = tramos(miss.values)
        longitudes = (g[:, 1] - g[:, 0] + 1) if len(g) else np.array([0])
        huecos.append(dict(estacion=k, minutos_faltantes=int(miss.sum()), pct=100 * miss.mean(), n_huecos=len(g),
                           mediana_min=float(np.median(longitudes)), p90_min=float(np.quantile(longitudes, .9)),
                           maximo_min=int(longitudes.max()),
                           hora_con_mas_huecos=int(miss.groupby(miss.index.hour).sum().idxmax()),
                           mes_con_mas_huecos=str(miss.groupby(miss.index.to_period("M")).sum().idxmax())))
        base = d.copy()
        v0 = ESTACIONES[k]["principal"]
        vv = base[v0].where(base[v0] != CENTINELA)
        gr = base.assign(_v=vv).groupby("calidad")
        t = gr.agg(registros=("calidad", "size"), primero=("calidad", lambda s: s.index.min()),
                   ultimo=("calidad", lambda s: s.index.max()))
        t["pct"] = 100 * t.registros / base["calidad"].notna().sum()
        rid = (base.calidad != base.calidad.shift()).cumsum()
        rach = base[base.calidad.notna()].groupby(rid).agg(c=("calidad", "first"), n=("calidad", "size"))
        t["rachas"] = rach.groupby("c").size()
        t["racha_mediana"] = rach.groupby("c").n.median()
        t["racha_max"] = rach.groupby("c").n.max()
        t["v_min"], t["v_mediana"], t["v_max"] = gr["_v"].min(), gr["_v"].median(), gr["_v"].max()
        t.insert(0, "variable_ref", v0)
        t.insert(0, "estacion", k)
        codigos.append(t.reset_index())
        s = serie_diaria(D, k).dropna()
        x = (s.index - s.index[0]).days.values
        lr, (tau, ptau) = stats.linregress(x, s.values), stats.kendalltau(x, s.values)
        kw = stats.kruskal(*[s[s.index.dayofweek == i].values for i in range(7)])
        tend.append(dict(estacion=k, variable=v0, dias=len(s), pendiente_por_30d=lr.slope * 30, p_regresion=lr.pvalue,
                         tau_kendall=tau, p_kendall=ptau, kruskal_dia_semana_p=kw.pvalue))
    pd.DataFrame(perfil).to_csv(os.path.join(out, "tabla_A1_perfil.csv"), index=False)
    pd.DataFrame(huecos).to_csv(os.path.join(out, "tabla_A2_huecos.csv"), index=False)
    pd.concat(codigos).to_csv(os.path.join(out, "tabla_A3_codigos_calidad.csv"), index=False)
    pd.DataFrame(tend).to_csv(os.path.join(out, "tabla_A4_tendencias.csv"), index=False)


def fig_estacion(D, k, out):
    cfg, c = ESTACIONES[k], COLOR[k]
    d, v = D.vals(k), ESTACIONES[k]["principal"]
    df, miss = D.df[k], D.faltantes(k)
    fig, ax = plt.subplots(2, 3, figsize=(11, 6.1))
    ax = ax.ravel()
    # A. serie diaria y tendencia
    s = serie_diaria(D, k)
    ax[0].plot(s.index, s.values, color=c, lw=.7, alpha=.55)
    ax[0].plot(s.rolling(7, min_periods=4).mean(), color=c, lw=1.8)
    sd = s.dropna()
    x = (sd.index - sd.index[0]).days.values
    lr = stats.linregress(x, sd.values)
    ax[0].plot(sd.index, lr.intercept + lr.slope * x, color="k", ls="--", lw=1)
    ax[0].set_title("A. Serie diaria y tendencia")
    ax[0].set_ylabel({"35": "Lluvia diaria p1 (mm)", "367": "Temperatura media (°C)", "803": "Nivel, mediana diaria",
                      "6004": "Radiación media (W/m²)"}[k])
    ax[0].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    # B. distribución
    z = d[v].dropna()
    b = ax[1]
    if k == "35":
        pulso = z[z > 0].min()
        vc = (z[z > 0] / pulso).round().astype(int).value_counts().sort_index()
        b.bar(vc.index * pulso, vc.values, width=pulso * .8, color=c)
        b.set_yscale("log")
        b.set_xlabel(f"p1 por minuto (mm); resolución {pulso:.3f}")
        b.set_title("B. Distribución de p1 (solo valores > 0)")
    elif k == "6004":
        b.hist(z[z > 0], bins=70, color=c)
        b.set_yscale("log")
        b.axvline(1361, color=ROJO, ls="--", lw=1)
        b.set_xlabel("Radiación (W/m²); línea roja: 1 361")
        b.set_title("B. Distribución de la radiación (> 0)")
    elif k == "803":
        b.hist(z, bins=80, color=c)
        b.set_yscale("log")
        b.axvline(0, color=ROJO, ls="--", lw=1)
        b.set_xlabel("Nivel (unidad por confirmar); línea roja: cero")
        b.set_title("B. Distribución del nivel")
    else:
        b.hist(z, bins=60, color=c)
        b.set_xlabel("Temperatura (°C)")
        b.set_title("B. Distribución de la temperatura")
    b.set_ylabel("Frecuencia (minutos)")
    # C. ciclo diurno
    cc = ax[2]
    if k == "35":
        h = d.p1.groupby(d.index.hour).sum()
        cc.bar(h.index, h.values, color=c)
        cc.set_ylabel("Lluvia acumulada (mm)")
        cc.set_title("C. Lluvia acumulada por hora del día")
    else:
        cd = ciclo_diurno(D, k)
        cc.fill_between(cd.index, cd.p25, cd.p75, color=c, alpha=.25)
        cc.plot(cd.index, cd.p50, color=c, lw=2, label="Marca de tiempo del archivo")
        if k == "6004":
            cs = ciclo_diurno(D, k, D.desfase_6004).reindex(range(24))
            cc.plot(cs.index, cs.p50, color="k", ls="--", lw=1.4, label=f"Corrida {D.desfase_6004:+d} h")
            cc.legend(frameon=False, fontsize=7)
        cc.set_ylabel({"367": "Temperatura (°C)", "803": "Nivel (mediana)", "6004": "Radiación (W/m²)"}[k])
        cc.set_title("C. Ciclo diurno (mediana y p25-p75)")
    cc.set_xlabel("Hora del día")
    cc.set_xticks(range(0, 24, 3))
    # D. huecos por hora
    h = miss.groupby(miss.index.hour).sum()
    ax[3].bar(h.index, h.values, color=GRIS)
    ax[3].set_title("D. Minutos faltantes por hora del día")
    ax[3].set_xlabel("Hora del día")
    ax[3].set_ylabel("Minutos")
    ax[3].set_xticks(range(0, 24, 3))
    # E. frecuencia de códigos
    ct = df.calidad.value_counts().sort_values()
    ax[4].barh([str(int(i)) for i in ct.index], ct.values, color=[GRIS if i == 1 else c for i in ct.index])
    ax[4].set_xscale("log")
    ax[4].set_title("E. Frecuencia de cada código")
    ax[4].set_xlabel("Registros (escala logarítmica)")
    for i, n in enumerate(ct.values):
        ax[4].text(n * 1.15, i, f"{n:,}".replace(",", " "), va="center", fontsize=7)
    ax[4].set_xlim(0.7, ct.max() * 8)
    # F. códigos distintos de 1 en el tiempo
    nn = df[df.calidad.notna() & (df.calidad != 1)]
    cods = sorted(nn.calidad.unique())
    pos = {c_: i for i, c_ in enumerate(cods)}
    ax[5].scatter(nn.index, nn.calidad.map(pos), s=3, color=c, alpha=.6, linewidths=0)
    ax[5].set_yticks(range(len(cods)))
    ax[5].set_yticklabels([str(int(c_)) for c_ in cods])
    ax[5].set_ylim(-.6, max(len(cods) - .4, .6))
    ax[5].set_title("F. Códigos distintos de 1 en el tiempo")
    ax[5].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    fig.suptitle(f"{cfg['nombre']}: resumen exploratorio", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, .96])
    fig.savefig(os.path.join(out, f"fig_A_{k}.png"))
    plt.close(fig)


def fig_367_distribuciones(D, out):
    d = D.vals("367")
    nm = {"h": "Humedad relativa (%)", "t": "Temperatura (°C)", "pr": "Presión (hPa)", "vv": "Velocidad del viento (m/s)",
          "vv_max": "Racha máxima (m/s)", "dv": "Dirección del viento (°)", "dv_max": "Dirección de la racha (°)",
          "p": "Precipitación por minuto (mm)"}
    fig, ax = plt.subplots(2, 4, figsize=(11, 4.6))
    for a, v in zip(ax.ravel(), ESTACIONES["367"]["vars"]):
        s = d[v].dropna()
        if v == "p":
            s = s[s > 0]
            a.set_yscale("log")
        a.hist(s, bins=60, color=COLOR["367"])
        a.set_xlabel(nm[v])
        a.set_title(v)
        if v == "h":
            a.set_yscale("log")
    ax[0, 0].set_ylabel("Frecuencia")
    ax[1, 0].set_ylabel("Frecuencia")
    fig.suptitle("Meteorológica 367: distribución de cada variable", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(os.path.join(out, "fig_A_367_distribuciones.png"))
    plt.close(fig)


def figs_comparativas(D, out):
    M = {k: D.faltantes(k) for k in ESTACIONES}
    nombres = [ESTACIONES[k]["nombre"] for k in ESTACIONES]
    # huecos por día y % de códigos distintos de 1
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6), gridspec_kw={"width_ratios": [1.7, 1]})
    dias = pd.date_range(D.ini, D.fin - pd.Timedelta(days=1), freq="D")
    mat = np.vstack([M[k].groupby(M[k].index.normalize()).sum().reindex(dias).fillna(0).values for k in ESTACIONES])
    im = ax[0].imshow(np.where(mat > 0, mat, np.nan), aspect="auto", cmap="YlOrRd", norm=LogNorm(1, max(mat.max(), 2)),
                      interpolation="nearest")
    ax[0].set_yticks(range(4))
    ax[0].set_yticklabels(nombres)
    tk = [i for i, d_ in enumerate(dias) if d_.day == 1]
    ax[0].set_xticks(tk)
    ax[0].set_xticklabels([dias[i].strftime("%b") for i in tk])
    ax[0].grid(False)
    ax[0].set_title("A. Minutos faltantes por día (escala logarítmica)")
    fig.colorbar(im, ax=ax[0], pad=.01, fraction=.03).set_label("min/día")
    pct = {k: 100 * D.etiqueta_debil(k).sum() / D.df[k]["calidad"].notna().sum() for k in ESTACIONES}
    ax[1].bar(list(ESTACIONES), [pct[k] for k in ESTACIONES], color=[COLOR[k] for k in ESTACIONES])
    for i, k in enumerate(ESTACIONES):
        ax[1].text(i, pct[k] + .3, f"{pct[k]:.1f} %", ha="center", fontsize=8)
    ax[1].set_title("B. Registros con calidad distinta de 1")
    ax[1].set_ylabel("% de registros")
    ax[1].set_xlabel("Estación")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_A_comp_huecos.png"))
    plt.close(fig)
    # ciclos diurnos normalizados
    def nrm(s):
        return (s - s.min()) / (s.max() - s.min())
    m, n, p, r = D.vals("367"), D.vals("803"), D.vals("35"), D.vals("6004")
    rad = r.radiacion.groupby(r.index.hour).mean().reindex(range(24)).fillna(0)
    rs = r.radiacion.copy()
    rs.index = rs.index + pd.Timedelta(hours=D.desfase_6004)
    rad_s = rs.groupby(rs.index.hour).mean().reindex(range(24)).fillna(0)
    fig, ax = plt.subplots(figsize=(8.2, 3.9))
    ax.plot(nrm(m.t.groupby(m.index.hour).mean()), color=COLOR["367"], lw=2, label="Temperatura 367")
    ax.plot(nrm(m.h.groupby(m.index.hour).mean()), color=COLOR["367"], lw=1.4, ls=":", label="Humedad 367")
    ax.plot(nrm(rad), color=COLOR["6004"], lw=2, label="Radiación 6004 (marca del archivo)")
    ax.plot(nrm(rad_s), color=COLOR["6004"], lw=1.6, ls="--", label=f"Radiación 6004 corrida {D.desfase_6004:+d} h")
    ax.plot(nrm(p.p1.groupby(p.index.hour).sum()), color=COLOR["35"], lw=1.6, label="Lluvia 35 (acumulado por hora)")
    ax.plot(nrm(n.nivel.groupby(n.index.hour).median()), color=COLOR["803"], lw=1.6, label="Nivel 803 (mediana)")
    ax.axvline(12, color=GRIS, lw=.8, ls="-.")
    ax.set_xticks(range(0, 24, 2))
    ax.set_xlabel("Hora del día según la marca de tiempo de cada archivo")
    ax.set_ylabel("Valor normalizado (0 a 1)")
    ax.set_title("Ciclos diurnos comparados")
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="upper left", bbox_to_anchor=(0, -.18))
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_A_comp_ciclos.png"))
    plt.close(fig)
    # lluvia 35 contra 367
    d35 = p.p1.resample("D").sum(min_count=1)
    d367 = m.p.resample("D").sum(min_count=1)
    ok = d35.notna() & d367.notna()
    lr = stats.linregress(d35[ok], d367[ok])
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    ax[0].scatter(d35[ok], d367[ok], s=12, color=COLOR["35"], alpha=.7)
    mx = max(d35.max(), d367.max())
    ax[0].plot([0, mx], [0, mx], color=GRIS, ls="--", lw=1, label="1 a 1")
    ax[0].plot([0, mx], [lr.intercept, lr.intercept + lr.slope * mx], color=ROJO, lw=1.5, label=f"Ajuste: pendiente {lr.slope:.2f}")
    ax[0].set_xlabel("Lluvia diaria en 35 (mm)")
    ax[0].set_ylabel("Lluvia diaria en 367 (mm)")
    ax[0].set_title("A. Lluvia diaria: estación 35 contra 367")
    ax[0].legend(frameon=False)
    ax[1].plot(d35.fillna(0).cumsum(), color=COLOR["35"], lw=2, label="35 (p1)")
    ax[1].plot(d367.fillna(0).cumsum(), color=COLOR["367"], lw=2, label="367 (p)")
    ax[1].set_title("B. Lluvia acumulada en el periodo")
    ax[1].set_ylabel("mm")
    ax[1].legend(frameon=False)
    ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_A_comp_lluvia.png"))
    plt.close(fig)
    # correlaciones horarias
    H = pd.DataFrame({
        "p1 (35)": p.p1.resample("h").sum(min_count=1), "p (367)": m.p.resample("h").sum(min_count=1),
        "t (367)": m.t.resample("h").mean(), "h (367)": m.h.resample("h").mean(),
        "pr (367)": m.pr.resample("h").mean(), "vv (367)": m.vv.resample("h").mean(),
        "nivel (803)": n.nivel.resample("h").mean(), "rad (6004) archivo": r.radiacion.resample("h").mean()})
    rh = r.radiacion.resample("h").mean()
    rh.index = rh.index + pd.Timedelta(hours=D.desfase_6004)
    H[f"rad (6004) {D.desfase_6004:+d} h"] = rh.reindex(H.index)
    C = H.corr(method="spearman")
    fig, ax = plt.subplots(figsize=(6.6, 5.4))
    im = ax.imshow(C.values, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.grid(False)
    ax.set_xticks(range(len(C)))
    ax.set_xticklabels(C.columns, rotation=40, ha="right", fontsize=7.5)
    ax.set_yticks(range(len(C)))
    ax.set_yticklabels(C.index, fontsize=7.5)
    for i in range(len(C)):
        for j in range(len(C)):
            ax.text(j, i, f"{C.values[i, j]:.2f}", ha="center", va="center", fontsize=6.5,
                    color="white" if abs(C.values[i, j]) > .55 else "black")
    fig.colorbar(im, fraction=.04, pad=.02)
    ax.set_title("Correlación de Spearman entre variables (datos horarios)")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_A_comp_correlacion.png"))
    plt.close(fig)
    # desfase horario del piranómetro
    fig, ax = plt.subplots(figsize=(6.2, 3))
    xs = sorted(D.corr_desfase)
    ax.plot(xs, [D.corr_desfase[s] for s in xs], marker="o", color=COLOR["6004"])
    ax.axvline(D.desfase_6004, color=ROJO, ls="--", lw=1)
    ax.set_xlabel("Desfase aplicado a la marca de tiempo (h)")
    ax.set_ylabel("Correlación con la geometría solar")
    ax.set_title(f"Desfase horario estimado del piranómetro: {D.desfase_6004:+d} h")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_A_desfase_piranometro.png"))
    plt.close(fig)


def parte_A(D, out):
    print("Tarea A: EDA")
    tablas_eda(D, out)
    for k in ESTACIONES:
        fig_estacion(D, k, out)
    fig_367_distribuciones(D, out)
    figs_comparativas(D, out)


# =====================================================================================
# 4. TAREA C1: REGLAS DETERMINISTAS
# =====================================================================================
def fuera(x, lo, hi):
    m = pd.Series(False, index=x.index)
    if lo is not None:
        m |= x < lo
    if hi is not None:
        m |= x > hi
    return m


def reglas_deterministas(D, k):
    """Devuelve un DataFrame booleano; cada columna es `tipo:nombre` y cada fila un minuto."""
    d, base, R = D.vals(k), D.df[k], {}
    R["F:faltante"] = ~base["presente"]
    R["F:centinela"] = D.sent[k]
    R["X:marca_irregular"] = base["seg_irregular"]
    R["D:duplicado"] = D.dup[k]
    # --- rango físico
    for v in ESTACIONES[k]["vars"]:
        lo, hi = LIMITES[k][v]
        if k == "6004":
            R["R:radiacion_sobre_cota"] = (d[v] > P["kt_max"] * D.g0 + 1.0) & (D.g0 > 0)
            R["R:radiacion_nocturna"] = (D.g0 <= 0) & (d[v] > P["noche_wm2"])
        R[f"R:{v}"] = fuera(d[v], lo, hi)
    # --- salto máximo
    for v, lim in SALTO_MAX.get(k, {}).items():
        R[f"J:{v}"] = d[v].diff().abs() > lim
    if k == "367":
        th = d.t.resample("h").mean()
        sh = th.diff().abs() >= SALTO_HORARIO_T
        R["J:t_horario"] = pd.Series(sh.reindex(d.index.floor("h")).fillna(False).astype(bool).values, index=d.index)
    # --- persistencia
    for v, n in CONGELADO_MIN.get(k, {}).items():
        x = d[v]
        if k == "6004":
            x = x.where((D.g0 > 50) & (x > 0))
        R[f"P:{v}"] = racha_mask(x, n)
    if k in ("35", "367"):
        for v in [c for c in ("p1", "p2", "p") if c in d.columns]:
            pulso = d[v][d[v] > 0].min()
            R[f"P:{v}_pulso_minimo"] = racha_mask(d[v].where(d[v] == pulso), P["lluvia_minima_min"])
    # --- coherencia interna
    if k == "367":
        R["I:vv_max_menor_vv"] = d.vv_max < d.vv
        R["I:calma_con_direccion"] = (d.vv == 0) & (d.dv > 0)
        R["I:lluvia_con_humedad_baja"] = (d.p > 0) & (d.h < 70)
        R["P:humedad100_sin_lluvia"] = racha_mask(d.h.where((d.h == 100) & (d.p == 0)), P["humedad100_min"])
    if k == "35":
        pulso = d.p1[d.p1 > 0].min()
        s1, s2 = d.p1.rolling(10, min_periods=10).sum(), d.p2.rolling(10, min_periods=10).sum()
        R["I:canales_discrepan_10min"] = (s1 - s2).abs() > P["discrep_pulsos"] * pulso
        d1, d2 = d.p1.resample("D").sum(min_count=1), d.p2.resample("D").sum(min_count=1)
        mal = ((d2 / d1.replace(0, np.nan) - 1).abs() > P["cociente_tol"]) & (np.maximum(d1, d2) >= P["cociente_min_mm"])
        mal_min = pd.Series(mal.reindex(d.index.floor("D")).fillna(False).astype(bool).values, index=d.index)
        R["I:cociente_diario_p2_p1"] = mal_min & ((d.p1 > 0) | (d.p2 > 0))
    out = pd.DataFrame(R).fillna(False).astype(bool)
    return out


def reglas_espaciales(D, R):
    """Reglas entre estaciones (tipo E)."""
    n = D.vals("803").nivel
    med = n.median()
    mad = 1.4826 * (n - med).abs().median()
    alto = (n - med) / mad > P["nivel_alto_z"]
    w = P["lluvia_previa_h"] * 60
    ll35 = D.vals("35").p1.fillna(0).rolling(w, min_periods=1).sum().shift(1)
    ll367 = D.vals("367").p.fillna(0).rolling(w, min_periods=1).sum().shift(1)
    R["803"]["E:nivel_alto_sin_lluvia_previa"] = alto & (ll35 < P["lluvia_previa_mm"]) & (ll367 < P["lluvia_previa_mm"])
    s35 = D.vals("35").p1.rolling(10, min_periods=10).sum()
    s367 = D.vals("367").p.rolling(10, min_periods=10).sum()
    R["35"]["E:lluvia_sin_respaldo_en_367"] = (s35 >= 2.0) & (s367 <= 0)
    return R


def familia(R, tipos="FXDRJPIE"):
    """Resume las columnas `tipo:nombre` en una columna por tipo."""
    return pd.DataFrame({t: R[[c for c in R.columns if c.startswith(t + ":")]].any(axis=1) for t in tipos
                         if any(c.startswith(t + ":") for c in R.columns)})


def reglas_sin_faltante(R):
    """Predicción de las reglas comparable con la etiqueta (excluye los minutos sin registro)."""
    cols = [c for c in R.columns if c != "F:faltante"]
    return R[cols].any(axis=1)


# =====================================================================================
# 5. TAREA C2: MÉTODOS ESTADÍSTICOS Y DE ML
# =====================================================================================
def descomposicion_robusta(s, ventana):
    """Tendencia (mediana móvil) + componente diurno (mediana por bloque de 10 min) + residuo."""
    tendencia = s.rolling(ventana, center=True, min_periods=max(10, ventana // 4)).median()
    det = s - tendencia
    bloque = s.index.hour * 6 + s.index.minute // 10
    estacional = det.groupby(bloque).transform("median")
    return tendencia, estacional, det - estacional


def z_robusto(res):
    """z robusto con mediana y MAD por hora del día; el MAD local no puede ser menor que el MAD global
    del residuo (evita que las horas muy estables o las series discretizadas disparen falsas alertas)."""
    hora = res.index.hour
    med = res.groupby(hora).transform("median")
    mad = (res - med).abs().groupby(hora).transform("median") * 1.4826
    mad_global = np.nanmedian((res - res.median()).abs()) * 1.4826
    mad = np.maximum(mad, mad_global).where(mad.notna())
    return (res - med) / mad.mask(mad <= 0)


def anomalias_mad(D, k):
    """Anomalías por residuos de la descomposición robusta. Devuelve (bool por variable, z, tendencia)."""
    d, flags, zs, tends = D.vals(k), {}, {}, {}
    if k == "35":
        return flags, zs, tends
    series = {}
    if k == "6004":
        series["radiacion"] = (d.radiacion / D.g0.where(D.g0 > 100))      # índice de claridad, solo de día
    else:
        for v in [c for c in ESTACIONES[k]["vars"] if c not in ("dv", "dv_max", "p")]:
            series[v] = d[v]
    for v, s in series.items():
        tend, _, res = descomposicion_robusta(s, P["mad_ventana"])
        z = z_robusto(res)
        flags[v], zs[v], tends[v] = (z.abs() > P["mad_k"]).fillna(False), z, tend
    return flags, zs, tends


def matriz_caracteristicas(D, k):
    d, idx = D.vals(k), D.idx
    hs = np.sin(2 * np.pi * (idx.hour + idx.minute / 60) / 24)
    hc = np.cos(2 * np.pi * (idx.hour + idx.minute / 60) / 24)
    F = {}
    if k == "35":
        F.update(p1=d.p1, p2=d.p2, dif=(d.p1 - d.p2).abs(), s10=d.p1.rolling(10, min_periods=1).sum(),
                 s60=d.p1.rolling(60, min_periods=1).sum())
    elif k == "367":
        for v in ["h", "t", "pr", "vv", "vv_max", "p"]:
            F[v] = d[v]
        for v in ["h", "t", "pr", "vv"]:
            F["d_" + v] = d[v].diff()
        for v in ["h", "t", "pr"]:
            F["sd10_" + v] = d[v].rolling(10, min_periods=5).std()
        F["dv_sin"], F["dv_cos"] = np.sin(np.deg2rad(d.dv)), np.cos(np.deg2rad(d.dv))
    elif k == "803":
        n = d.nivel
        F.update(nivel=n, d1=n.diff(), d10=n.diff(10), sd10=n.rolling(10, min_periods=5).std(),
                 dev60=n - n.rolling(61, center=True, min_periods=20).median())
    else:
        g = D.g0.where(D.g0 > 50)
        kt = d.radiacion / g
        F.update(kt=kt, d_kt=kt.diff(), sd10=kt.rolling(10, min_periods=5).std(), rad=d.radiacion.where(g.notna()))
    F["hs"], F["hc"] = pd.Series(hs, index=idx), pd.Series(hc, index=idx)
    return pd.DataFrame(F)


def isolation_forest(D, k):
    X = matriz_caracteristicas(D, k).replace([np.inf, -np.inf], np.nan).dropna()
    Xs = RobustScaler().fit_transform(X.values)
    modelo = IsolationForest(n_estimators=P["if_arboles"], contamination=P["if_contaminacion"],
                             max_samples=min(P["if_max_muestras"], len(X)), random_state=P["if_semilla"], n_jobs=-1)
    modelo.fit(Xs)
    score = pd.Series(modelo.score_samples(Xs), index=X.index).reindex(D.idx)
    flag = pd.Series(modelo.predict(Xs) == -1, index=X.index).reindex(D.idx).fillna(False).astype(bool)
    return flag, score


def cambios_de_regimen(s, ventana_h, k):
    """Desplazamiento de la mediana horaria (tras quitar el ciclo diurno) entre ventanas contiguas."""
    h = s.resample("h").median()
    h = h - h.groupby(h.index.hour).transform("median")
    izq = h.rolling(ventana_h, min_periods=ventana_h // 2).median().shift(1)
    der = h[::-1].rolling(ventana_h, min_periods=ventana_h // 2).median()[::-1]
    delta = der - izq
    sc = 1.4826 * np.nanmedian(np.abs(delta - np.nanmedian(delta)))
    if not np.isfinite(sc) or sc == 0:
        sc = np.nanstd(delta)
    z = (delta.abs() / sc).fillna(0)
    pk, _ = find_peaks(z.values, height=k, distance=ventana_h)
    return pd.DataFrame({"t": h.index[pk], "delta": delta.values[pk], "z": z.values[pk]})


def regimenes(D, k):
    series = {"367": ["t", "h", "pr"], "803": ["nivel"], "6004": ["radiacion"]}.get(k, [])
    d = D.vals(k)
    mascara = pd.Series(False, index=D.idx)
    cps = []
    for v in series:
        c = cambios_de_regimen(d[v], P["regimen_ventana_h"], P["regimen_k"])
        c.insert(0, "variable", v)
        c.insert(0, "estacion", k)
        cps.append(c)
        for t in c.t:
            mascara.loc[t - pd.Timedelta(minutes=30): t + pd.Timedelta(minutes=30)] = True
    return mascara, (pd.concat(cps) if cps else pd.DataFrame(columns=["estacion", "variable", "t", "delta", "z"]))


# =====================================================================================
# 6. TAREA C3: COMPARACIÓN
# =====================================================================================
def metricas(pred, lab, evaluable, tol):
    p, l = pred[evaluable].values.astype(bool), lab[evaluable].values.astype(bool)
    tp = (p & l).sum()
    prec, rec = tp / max(p.sum(), 1), tp / max(l.sum(), 1)
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    dp, dl = dilatar(p, tol), dilatar(l, tol)
    prec_t, rec_t = (p & dl).sum() / max(p.sum(), 1), (l & dp).sum() / max(l.sum(), 1)
    f1_t = 2 * prec_t * rec_t / (prec_t + rec_t) if prec_t + rec_t else 0.0
    return dict(n_alertas=int(p.sum()), n_etiqueta=int(l.sum()), precision=prec, exhaustividad=rec, f1=f1,
                precision_tol=prec_t, exhaustividad_tol=rec_t, f1_tol=f1_t)


def comparar(D, M):
    """M[k] = dict de máscaras: reglas, mad, iforest, regimen, estadistico, union, etiqueta."""
    filas, acuerdo, por_codigo, desacuerdos = [], [], [], []
    for k in ESTACIONES:
        m, ev = M[k], D.df[k]["calidad"].notna()
        for nombre in ["reglas", "mad", "iforest", "regimen", "estadistico", "union"]:
            filas.append(dict(estacion=k, metodo=nombre, **metricas(m[nombre], m["etiqueta"], ev, P["tolerancia_min"])))
        # acuerdo entre enfoques
        a, b = m["reglas"][ev].values, m["estadistico"][ev].values
        acuerdo.append(dict(estacion=k, solo_reglas=int((a & ~b).sum()), solo_estadistico=int((~a & b).sum()),
                            ambos=int((a & b).sum()), ninguno=int((~a & ~b).sum()),
                            jaccard=(a & b).sum() / max((a | b).sum(), 1),
                            kappa=cohen_kappa_score(a, b) if (a.any() and b.any()) else np.nan))
        # exhaustividad por código de calidad
        cal = D.df[k]["calidad"]
        for c in sorted(cal.dropna().unique()):
            s = cal == c
            fila = dict(estacion=k, codigo=int(c), registros=int(s.sum()))
            for nombre in ["reglas", "mad", "iforest", "regimen", "estadistico", "union"]:
                fila[nombre] = float(m[nombre][s].mean())
            por_codigo.append(fila)
        # casos de desacuerdo agrupados en eventos
        L, R_, S = m["etiqueta"] & ev, m["reglas"] & ev, m["estadistico"] & ev
        categorias = {"solo_etiqueta": L & ~R_ & ~S, "solo_reglas": R_ & ~S & ~L, "solo_estadistico": S & ~R_ & ~L,
                      "reglas_y_estadistico_sin_etiqueta": R_ & S & ~L}
        for cat, mk in categorias.items():
            for a_, b_ in tramos(mk.values, P["brecha_eventos"]):
                seg = D.idx[a_:b_ + 1]
                desacuerdos.append(dict(estacion=k, categoria=cat, inicio=seg[0], fin=seg[-1], minutos=len(seg),
                                        codigos=",".join(str(int(c)) for c in sorted(cal.iloc[a_:b_ + 1].dropna().unique())),
                                        reglas=int(R_.iloc[a_:b_ + 1].sum()), estadistico=int(S.iloc[a_:b_ + 1].sum())))
    return pd.DataFrame(filas), pd.DataFrame(acuerdo), pd.DataFrame(por_codigo), pd.DataFrame(desacuerdos)


# =====================================================================================
# 7. TAREA C4: CLASIFICACIÓN DE ANOMALÍAS
# =====================================================================================
CLASES = ["falla_sensor", "falla_comunicacion", "configuracion", "evento_fisico_real", "indeterminado"]
ACCION = {"falla_sensor": "excluir de productos y programar revisión del sensor",
          "falla_comunicacion": "no imputar a ciegas; recuperar del datalogger",
          "configuracion": "corregir metadatos o configuración y reprocesar",
          "evento_fisico_real": "conservar, marcar y no eliminar",
          "indeterminado": "conservar marcado hasta contrastar con otras fuentes"}


def hallazgos_configuracion(D):
    """Problemas de configuración que se ven en toda la serie y no en un evento puntual."""
    H = []
    mejor, corr = D.desfase_6004, D.corr_desfase
    if mejor != 0:
        H.append(dict(estacion="6004", hallazgo="zona_horaria",
                      evidencia=f"correlación con la geometría solar {corr[mejor]:.3f} con desfase {mejor:+d} h frente a {corr[0]:.3f} sin corregir",
                      recomendacion="declarar la zona horaria; guardar en UTC y convertir en los productos"))
    p, m = D.vals("35").p1.resample("D").sum(min_count=1), D.vals("367").p.resample("D").sum(min_count=1)
    ok = p.notna() & m.notna() & ((p > 0) | (m > 0))
    if ok.sum() > 10:
        lr = stats.linregress(p[ok], m[ok])
        if abs(lr.rvalue) > .95 and abs(lr.slope - 1) > .2:
            H.append(dict(estacion="35/367", hallazgo="escala_de_lluvia",
                          evidencia=f"r={lr.rvalue:.3f} y pendiente {lr.slope:.2f} entre acumulados diarios (35 contra 367)",
                          recomendacion="revisar factor de escala, resolución y captación de uno de los pluviómetros"))
    n = D.vals("803").nivel.dropna()
    if (n < 0).mean() > .01:
        H.append(dict(estacion="803", hallazgo="nivel_negativo_cota_cero",
                      evidencia=f"{(n < 0).mean() * 100:.1f} % de valores negativos (mínimo {n.min():.1f})",
                      recomendacion="documentar unidad y cota cero de la mira; definir el rango por estación"))
    d367 = D.vals("367")
    if ((d367.dv == 360).sum() == 0) and ((d367.dv_max == 360).sum() > 0):
        H.append(dict(estacion="367", hallazgo="codificacion_del_norte",
                      evidencia=f"dv no usa 360 y dv_max sí ({int((d367.dv_max == 360).sum())} casos); ceros: dv {int((d367.dv == 0).sum())}, dv_max {int((d367.dv_max == 0).sum())}",
                      recomendacion="normalizar la convención de la dirección en la capa validada"))
    paso35 = D.vals("35").p1[D.vals("35").p1 > 0].min()
    paso367 = d367.p[d367.p > 0].min()
    if abs(paso35 - paso367) / max(paso35, paso367) > .5:
        H.append(dict(estacion="35/367", hallazgo="resolucion_de_lluvia_distinta",
                      evidencia=f"resolución mínima {paso35:.3f} mm en 35 y {paso367:.3f} mm en 367",
                      recomendacion="declarar el tipo de sensor y la resolución en los metadatos"))
    for k in ESTACIONES:
        n_irr = int(D.df[k]["seg_irregular"].sum())
        if n_irr:
            H.append(dict(estacion=k, hallazgo="marcas_de_tiempo_con_segundos", evidencia=f"{n_irr} registros con segundos distintos de cero",
                          recomendacion="alinear las marcas al minuto en la ingesta"))
    return pd.DataFrame(H)


def clasificar_eventos(D, M, R, Z, TEND):
    """Agrupa los minutos relevantes en eventos y asigna una clase con su justificación."""
    ll35 = D.vals("35").p1.fillna(0).rolling(P["lluvia_previa_h"] * 60, min_periods=1).sum()
    ll367 = D.vals("367").p.fillna(0).rolling(P["lluvia_previa_h"] * 60, min_periods=1).sum()
    com = {k: dilatar(D.faltantes(k).values, 2) for k in ESTACIONES}
    acum = {k: np.r_[0, np.cumsum(com[k])] for k in com}
    filas = []
    for k in ESTACIONES:
        m, F = M[k], familia(R[k])
        rel = (m["reglas"] & (F.get("F", False) | F.get("R", False) | F.get("P", False) | F.get("X", False))) | \
              ((m["reglas"].astype(int) + m["estadistico"].astype(int) + m["etiqueta"].astype(int)) >= 2)
        rel = rel | R[k]["F:faltante"]
        d, cal, vdom = D.vals(k), D.df[k]["calidad"], ESTACIONES[k]["principal"]
        for a, b in tramos(rel.values, P["brecha_eventos"]):
            sl = slice(a, b + 1)
            dur = b - a + 1
            tags = {t: bool(F[t].iloc[sl].any()) for t in F.columns}
            cols = R[k].iloc[sl].sum()
            cols = cols[cols > 0]
            regla_dom = cols.idxmax() if len(cols) else ""
            var = regla_dom.split(":")[1] if ":" in regla_dom and regla_dom.split(":")[1] in d.columns else vdom
            x = d[var].iloc[sl]
            ref = TEND.get(k, {}).get(var)
            if ref is not None:
                dev = (x - ref.iloc[sl]).abs()
                pico = float(dev.max()) if dev.notna().any() else np.nan
            else:
                pico = float(x.abs().max()) if x.notna().any() else np.nan
            cods = cal.iloc[sl].dropna()
            cod = int(cods.mode().iloc[0]) if len(cods) and (cods != 1).any() else (int(cods.iloc[0]) if len(cods) else -1)
            if (cods != 1).any():
                cod = int(cods[cods != 1].mode().iloc[0])
            otras = sum(int(acum[o][min(b + 3, len(D.idx))] - acum[o][max(a - 2, 0)] > 0) for o in ESTACIONES if o != k)
            lluvia = float(max(ll35.iloc[a], ll367.iloc[a]))
            fuente = ",".join(n for n, c in [("reglas", m["reglas"]), ("mad", m["mad"]), ("iforest", m["iforest"]),
                                             ("regimen", m["regimen"]), ("etiqueta", m["etiqueta"])] if c.iloc[sl].any())
            clase, conf, why = "indeterminado", "baja", f"señalado por {fuente or 'ninguna fuente'} sin una regla determinista que identifique la causa"
            if tags.get("F") or tags.get("X"):
                clase = "falla_comunicacion"
                if otras:
                    conf, why = "alta", f"faltante o centinela que coincide con huecos en {otras} estación(es) más: causa común en la recepción"
                else:
                    conf, why = "media", "faltante o centinela sin coincidencia en otras estaciones"
            elif tags.get("P"):
                clase, conf, why = "falla_sensor", "alta", f"valor congelado ({regla_dom}) durante {dur} min"
            elif k == "6004" and tags.get("R") and "nocturna" in regla_dom:
                clase, conf, why = "configuracion", "alta", "radiación en horas sin sol: zona horaria o reloj mal configurado"
            elif k == "6004" and tags.get("R"):
                kt_pico = float((d.radiacion.iloc[sl] / D.g0.iloc[sl].where(D.g0.iloc[sl] > 0)).max())
                if kt_pico <= 1.3 and dur >= 2:
                    clase, conf, why = "evento_fisico_real", "media", f"índice de claridad máximo {kt_pico:.2f} con persistencia: posible realce por nubes"
                else:
                    clase, conf, why = "falla_sensor", "media", f"índice de claridad {kt_pico:.2f} en {dur} min, sin respaldo físico"
            elif k == "803" and tags.get("R") and (x < 0).mean() > .5:
                clase, conf, why = "configuracion", "media", "nivel negativo sostenido: posible cota cero o desplazamiento del cero"
            elif k == "803" and (tags.get("E") or (x.max() > 10 * max(abs(d.nivel.median()), 1))):
                if lluvia >= P["lluvia_previa_mm"] and dur >= 30:
                    clase, conf, why = "evento_fisico_real", "media", f"nivel alto sostenido {dur} min con {lluvia:.1f} mm de lluvia previa"
                elif dur >= 30:
                    clase, conf, why = "indeterminado", "baja", f"nivel alto sostenido {dur} min sin lluvia previa en 35 ni 367: requiere estación aguas arriba"
                else:
                    clase, conf, why = "falla_sensor", "media", f"pico aislado de {dur} min sin respaldo"
            elif tags.get("J") and dur <= 3:
                clase, conf, why = "falla_sensor", "media", f"salto aislado ({regla_dom}) de {dur} min que revierte"
            elif tags.get("J") and k == "367" and lluvia > 0 and float(d.p.iloc[max(a - 30, 0):b + 30].sum()) > 0:
                clase, conf, why = "evento_fisico_real", "media", "salto de temperatura acompañado de lluvia en la misma ventana"
            elif tags.get("I") and k == "35":
                p1, p2 = d.p1.iloc[sl].sum(), d.p2.iloc[sl].sum()
                if min(p1, p2) > 0:
                    clase, conf, why = "evento_fisico_real", "media", "ambos canales registran lluvia; la diferencia es de intensidad"
                else:
                    clase, conf, why = "falla_sensor", "baja", "un canal registra lluvia y el otro no"
            elif tags.get("I") and k == "367":
                clase, conf, why = "falla_sensor", "baja", f"incoherencia entre variables ({regla_dom})"
            elif k == "35" and dur >= 500 and m["etiqueta"].iloc[sl].mean() > .9:
                clase, conf, why = "configuracion", "baja", f"tramo de {dur} min marcado en bloque por el equipo con valores sin anomalía: posible mantenimiento"
            elif k == "35" and m["iforest"].iloc[sl].any() and (d.p1.iloc[sl] > 0).any() and (d.p2.iloc[sl] > 0).any():
                clase, conf, why = "evento_fisico_real", "media", "intensidad extrema registrada por los dos canales"
            elif m["regimen"].iloc[sl].any() and not m["reglas"].iloc[sl].any():
                clase, conf, why = "configuracion", "baja", "cambio persistente de nivel sin corroboración: posible ajuste o calibración"
            filas.append(dict(estacion=k, variable=var, inicio=D.idx[a], fin=D.idx[b], minutos=dur, fuentes=fuente,
                              reglas=";".join(cols.index[:3]), codigo_equipo=cod, valor_pico_desviacion=pico,
                              lluvia_previa_mm=round(lluvia, 2), coincide_en_otras=otras, clase=clase,
                              confianza=conf, justificacion=why, accion=ACCION[clase],
                              proteger_extremo=(clase == "evento_fisico_real")))
    return pd.DataFrame(filas)


# =====================================================================================
# 8. TAREA C5: BANDERAS Y DICCIONARIO DE CÓDIGOS
# =====================================================================================
BANDERAS = pd.DataFrame([
    # codigo, nombre, familia, origen, accion
    ("OK", "Aprobado", "ok", "todas las pruebas", "entra a los productos"),
    ("F", "Faltante", "rechazo", "regla", "no entra a agregados; cuenta contra la cobertura"),
    ("X", "Formato inválido", "rechazo", "regla", "se excluye; se conserva el crudo"),
    ("D", "Duplicado", "rechazo", "regla", "se conserva uno; el resto va a cuarentena"),
    ("R", "Fuera de rango", "rechazo", "regla", "se excluye con causa registrada"),
    ("J", "Salto máximo", "sospechoso", "regla", "se conserva marcado para revisión"),
    ("P", "Persistencia o valor congelado", "sospechoso", "regla", "se conserva marcado; excluir si se confirma falla"),
    ("A", "Anomalía estadística (MAD)", "sospechoso", "estadístico", "se conserva marcado para revisión"),
    ("N", "Anomalía multivariada (Isolation Forest)", "sospechoso", "ML", "se conserva marcado para revisión"),
    ("T", "Cambio de régimen", "sospechoso", "estadístico", "se conserva marcado; revisar cota cero o calibración"),
    ("I", "Incoherente", "incoherente", "regla", "se conserva marcado para revisión humana"),
    ("E", "Inconsistencia espacial", "incoherente", "regla", "revisión humana, sin rechazo automático"),
    ("Q", "Marcado por el equipo (calidad distinta de 1)", "contexto", "etiqueta débil", "informativo; no decide por sí solo"),
    ("M", "Imputado", "modificado", "tratamiento", "se publica solo si el producto lo permite, siempre marcado"),
    ("C", "Corregido", "modificado", "revisión", "se conserva el original junto al corregido"),
    ("K", "Mantenimiento", "contexto", "bitácora", "se excluye de productos y de la cobertura"),
], columns=["codigo", "bandera", "familia", "origen", "que_pasa_con_el_dato"])
BANDERAS["bit"] = [0 if c == "OK" else i for i, c in enumerate(BANDERAS.codigo)]
PRIORIDAD = ["F", "X", "D", "R", "E", "I", "J", "P", "A", "N", "T", "Q"]


def asignar_banderas(D, k, R, M, Z):
    """Una columna booleana por bandera; devuelve el DataFrame de banderas, la cadena y la principal."""
    F = familia(R)
    B = pd.DataFrame(False, index=D.idx, columns=[c for c in BANDERAS.codigo if c != "OK"])
    for t in "FXDRJPIE":
        if t in F:
            B[t] = F[t]
    B["A"] = M["mad"]
    B["N"] = M["iforest"]
    B["T"] = M["regimen"]
    B["Q"] = M["etiqueta"]
    cadena = B.apply(lambda c: np.where(c, c.name, ""), axis=0)
    texto = cadena.apply(lambda f: ",".join([x for x in f if x]) or "OK", axis=1)
    principal = pd.Series("OK", index=D.idx)
    for c in reversed(PRIORIDAD):
        principal = principal.mask(B[c], c)
    peso = BANDERAS.set_index("codigo")["bit"]
    bits = sum(B[c].astype(np.int64) * (1 << int(peso[c])) for c in B.columns)
    return B, texto, principal, bits


def inferir_codigos(D, k, R, M):
    """Diccionario de códigos de la columna `calidad` a partir de las firmas que dejan en los datos."""
    d, cal = D.vals(k), D.df[k]["calidad"]
    raw = D.df[k]
    F = familia(R)
    v0 = ESTACIONES[k]["principal"]
    base = cal == 1
    mediana1 = d[v0][base].median()
    p99_1 = d[v0][base].quantile(.99)
    filas = []
    rid = (cal != cal.shift()).cumsum()
    rachas = cal[cal.notna()].groupby(rid).agg(["first", "size"])
    for c in sorted(cal.dropna().unique()):
        s = cal == c
        r = rachas[rachas["first"] == c]["size"]
        x = d[v0][s]
        firmas = []
        if c == 1:
            filas.append(dict(estacion=k, codigo=1, registros=int(s.sum()), pct=100 * s.sum() / cal.notna().sum(),
                              rachas=len(r), racha_mediana=float(r.median()), racha_max=int(r.max()),
                              firmas="mayoría de los registros", significado_inferido="dato normal (suposición)",
                              bandera_equivalente="OK", confianza="media"))
            continue
        if k == "35":
            sent = (raw.loc[s, ["p1", "p2"]] == CENTINELA).all(axis=1).mean()
            if sent > .5:
                firmas.append(f"{sent:.0%} de los registros trae el centinela -999")
            c1 = ((raw.p1 > 0) & (raw.p2 == 0))[s].mean()
            c2 = ((raw.p2 > 0) & (raw.p1 == 0))[s].mean()
            if c1 > .9:
                firmas.append("p1 > 0 con p2 = 0 (el canal 2 no registra)")
            if c2 > .9:
                firmas.append("p2 > 0 con p1 = 0 (el canal 1 no registra)")
        if k == "367":
            if (d.h[s] == 100).mean() > .9:
                firmas.append("humedad en 100 % en todo el tramo")
            pos = d.p[s]
            if pos.nunique() == 1 and pos.iloc[0] > 0:
                firmas.append(f"lluvia constante de {pos.iloc[0]:.2f} mm en todos los registros")
        if k == "803":
            if (x < 0).mean() > .9:
                firmas.append("valores negativos en casi todos los registros")
            if x.median() > 10 * max(abs(mediana1), 1):
                firmas.append(f"mediana {x.median():.1f}, más de 10 veces la del código 1")
        if k == "6004":
            g = D.g0[s]
            if (x > 1361).any() or x.median() > p99_1 * .8:
                firmas.append("valores muy altos (picos) con cambios bruscos entre minutos")
            if ((g > 0) & (g < 150)).mean() > .5 and x.median() < .05 * max(p99_1, 1):
                firmas.append("valores muy bajos en horas de amanecer o atardecer")
        if len(r) and r.max() >= 500 and r.max() >= .8 * s.sum():
            firmas.append(f"un solo tramo continuo de {int(r.max())} min (marcado en bloque)")
        elif len(r) and r.median() <= 2:
            firmas.append("rachas cortas y dispersas")
        for nombre, f_ in (("regla", M["reglas"]), ("estadístico", M["estadistico"])):
            pct = float(f_[s].mean())
            if pct >= .5:
                firmas.append(f"{pct:.0%} de sus registros lo detecta también el enfoque {nombre}")
        sig, bandera, conf = "sin patrón claro en las columnas", "Q", "baja"
        txt = " ".join(firmas)
        if "centinela" in txt:
            sig, bandera, conf = "dato faltante (centinela -999)", "F", "alta"
        elif "canal 2 no registra" in txt:
            sig, bandera, conf = "discrepancia entre canales: solo registra el canal 1", "I", "alta"
        elif "canal 1 no registra" in txt:
            sig, bandera, conf = "discrepancia entre canales: solo registra el canal 2", "I", "alta"
        elif "humedad en 100" in txt:
            sig, bandera, conf = "humedad congelada en saturación", "P", "media"
        elif "lluvia constante" in txt:
            sig, bandera, conf = "lluvia mínima aislada (un solo paso)", "J", "media"
        elif "valores negativos" in txt:
            sig, bandera, conf = "nivel negativo (cota cero o desplazamiento)", "R", "media"
        elif "más de 10 veces" in txt:
            sig, bandera, conf = "nivel alto o extremo marcado por el equipo", "J", "media"
        elif "picos" in txt:
            sig, bandera, conf = "picos de radiación (posible exceso sobre la cota)", "R", "media"
        elif "amanecer o atardecer" in txt:
            sig, bandera, conf = "radiación casi nula en horas de transición", "J", "baja"
        elif "marcado en bloque" in txt:
            sig, bandera, conf = "periodo marcado en bloque (mantenimiento o fallo prolongado)", "K", "baja"
        filas.append(dict(estacion=k, codigo=int(c), registros=int(s.sum()), pct=100 * s.sum() / cal.notna().sum(),
                          rachas=len(r), racha_mediana=float(r.median()), racha_max=int(r.max()),
                          firmas="; ".join(firmas) if firmas else "ninguna", significado_inferido=sig,
                          bandera_equivalente=bandera, confianza=conf))
    return pd.DataFrame(filas)


# =====================================================================================
# 9. GRÁFICAS DE LA TAREA C
# =====================================================================================
def fig_C1(D, R, out):
    fig, ax = plt.subplots(2, 2, figsize=(11, 6.4))
    for a, k in zip(ax.ravel(), ESTACIONES):
        cnt = R[k].sum().sort_values()
        cnt = cnt[cnt > 0]
        a.barh(cnt.index, cnt.values, color=COLOR[k])
        a.set_xscale("log")
        a.set_title(f"{ESTACIONES[k]['nombre']}: minutos señalados por regla")
        a.set_xlabel("Minutos (escala logarítmica)")
        a.tick_params(axis="y", labelsize=7)
        for i, n in enumerate(cnt.values):
            a.text(n * 1.15, i, f"{n:,}".replace(",", " "), va="center", fontsize=7)
        a.set_xlim(0.7, max(cnt.max() * 10, 10) if len(cnt) else 10)
    fig.suptitle("C1. Reglas deterministas", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, .96])
    fig.savefig(os.path.join(out, "fig_C1_reglas.png"))
    plt.close(fig)


def fig_C2(D, M, Z, cps, out):
    fig, ax = plt.subplots(4, 1, figsize=(11, 9.4))
    for a, k in zip(ax, ESTACIONES):
        v = ESTACIONES[k]["principal"]
        s = D.vals(k)[v]
        h = s.resample("h").sum(min_count=1) if k == "35" else s.resample("h").mean()
        a.plot(h.index, h.values, color=COLOR[k], lw=.8, alpha=.8, label=f"{v} (horario)")
        for nombre, marca, col in [("mad", "o", "#d69e2e"), ("iforest", "x", ROJO), ("regimen", "|", "k")]:
            fl = M[k][nombre]
            if fl.any():
                hh = fl.astype(int).resample("h").max().astype(bool)
                horas = hh[hh].index
                a.scatter(horas, h.reindex(horas), marker=marca, s=12 if marca != "|" else 70, color=col, linewidths=.8,
                          alpha=.55 if nombre == "mad" else .9,
                          label={"mad": "Residuo MAD", "iforest": "Isolation Forest", "regimen": "Cambio de régimen"}[nombre])
        a.set_title(f"{ESTACIONES[k]['nombre']}")
        a.legend(frameon=False, fontsize=7, ncol=4, loc="upper right")
        a.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    fig.suptitle("C2. Anomalías estadísticas y de aprendizaje automático sobre la variable principal", x=0.01, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, .97])
    fig.savefig(os.path.join(out, "fig_C2_metodos.png"))
    plt.close(fig)


def fig_C3(metricas_df, acuerdo_df, codigos_df, out):
    metodos = ["reglas", "mad", "iforest", "regimen", "estadistico", "union"]
    cols = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#d69e2e", "#6b7280"]
    fig, ax = plt.subplots(2, 2, figsize=(11, 6.4))
    for a, (m_, t) in zip(ax.ravel(), [("precision", "Precisión"), ("exhaustividad", "Exhaustividad"),
                                       ("precision_tol", "Precisión con tolerancia"), ("exhaustividad_tol", "Exhaustividad con tolerancia")]):
        piv = metricas_df.pivot(index="estacion", columns="metodo", values=m_).reindex(list(ESTACIONES))[metodos]
        w = .13
        for i, mt in enumerate(metodos):
            a.bar(np.arange(len(piv)) + (i - 2.5) * w, piv[mt].values, width=w, color=cols[i], label=mt)
        a.set_xticks(range(len(piv)))
        a.set_xticklabels(piv.index)
        a.set_ylim(0, 1.05)
        a.set_title(t + " frente a la columna calidad")
    ax[0, 0].legend(frameon=False, fontsize=7, ncol=3)
    fig.suptitle(f"C3. Métodos frente a la etiqueta débil (tolerancia ±{P['tolerancia_min']} min)", x=0.01, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, .96])
    fig.savefig(os.path.join(out, "fig_C3_metricas.png"))
    plt.close(fig)
    # exhaustividad por código
    cod = codigos_df.copy()
    cod["fila"] = cod.estacion + ": " + cod.codigo.astype(str)
    mat = cod.set_index("fila")[metodos]
    fig, ax = plt.subplots(figsize=(6.6, max(3.5, .28 * len(mat) + 1.2)))
    im = ax.imshow(mat.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.grid(False)
    ax.set_xticks(range(len(metodos)))
    ax.set_xticklabels(metodos, rotation=30, ha="right")
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels(mat.index, fontsize=7)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat.values[i, j]:.2f}", ha="center", va="center", fontsize=6.5,
                    color="white" if mat.values[i, j] > .6 else "black")
    fig.colorbar(im, fraction=.04, pad=.02).set_label("Proporción de registros detectados")
    ax.set_title("Detección por código de calidad (estación: código)")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_C3_codigos.png"))
    plt.close(fig)
    # acuerdo entre enfoques
    fig, ax = plt.subplots(figsize=(7, 3.6))
    a_ = acuerdo_df.set_index("estacion").reindex(list(ESTACIONES))
    bottom = np.zeros(len(a_))
    for col, lab, cc in [("ambos", "Ambos", "#1baf7a"), ("solo_reglas", "Solo reglas", "#2a78d6"),
                         ("solo_estadistico", "Solo estadístico", "#eb6834")]:
        ax.bar(a_.index, a_[col].values, bottom=bottom, label=lab, color=cc)
        bottom += a_[col].values
    ax.set_yscale("symlog")
    ax.set_ylabel("Minutos señalados")
    ax.set_title("Acuerdo entre reglas y métodos estadísticos")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_C3_acuerdo.png"))
    plt.close(fig)


def fig_C4(eventos, out):
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.8), gridspec_kw={"width_ratios": [1, 1, 1.5]})
    cc = {"falla_sensor": "#e34948", "falla_comunicacion": "#d69e2e", "configuracion": "#6250d6",
          "evento_fisico_real": "#1baf7a", "indeterminado": "#9ca3af"}
    cnt = eventos.groupby(["estacion", "clase"]).size().unstack(fill_value=0).reindex(list(ESTACIONES)).reindex(columns=CLASES).fillna(0)
    mins = eventos.groupby(["estacion", "clase"]).minutos.sum().unstack(fill_value=0).reindex(list(ESTACIONES)).reindex(columns=CLASES).fillna(0)
    for a, tabla, tit in [(ax[0], cnt, "Eventos por clase"), (ax[1], mins, "Minutos por clase")]:
        b = np.zeros(len(tabla))
        for c in CLASES:
            a.bar(tabla.index, tabla[c].values, bottom=b, color=cc[c], label=c)
            b += tabla[c].values
        a.set_title(tit)
        a.set_xlabel("Estación")
    ax[1].set_yscale("log")
    ax[0].legend(frameon=False, fontsize=6.5)
    y = {k: i for i, k in enumerate(ESTACIONES)}
    ax[2].scatter(eventos.inicio, eventos.estacion.map(y), s=np.clip(eventos.minutos, 3, 400) / 4 + 4,
                  c=eventos.clase.map(cc), alpha=.6, linewidths=0)
    ax[2].set_yticks(list(y.values()))
    ax[2].set_yticklabels(list(y))
    ax[2].set_title("Eventos en el tiempo (tamaño: duración)")
    ax[2].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    fig.suptitle("C4. Clasificación de anomalías relevantes", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(os.path.join(out, "fig_C4_clases.png"))
    plt.close(fig)


def fig_C5(B, out):
    fig, ax = plt.subplots(figsize=(10, 3.8))
    codigos = [c for c in BANDERAS.codigo if c != "OK"]
    w = .2
    for i, k in enumerate(ESTACIONES):
        ax.bar(np.arange(len(codigos)) + (i - 1.5) * w, [int(B[k][c].sum()) if c in B[k] else 0 for c in codigos], width=w,
               color=COLOR[k], label=k)
    ax.set_yscale("symlog")
    ax.set_xticks(range(len(codigos)))
    ax.set_xticklabels(codigos)
    ax.set_ylabel("Minutos con la bandera (escala symlog)")
    ax.set_title("C5. Banderas asignadas por estación")
    ax.legend(frameon=False, ncol=4)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_C5_banderas.png"))
    plt.close(fig)


# =====================================================================================
# 10. ORQUESTACIÓN DE LA TAREA C
# =====================================================================================
def parte_C(D, out, exportar):
    print("Tarea C: reglas (C1)")
    R = {k: reglas_deterministas(D, k) for k in ESTACIONES}
    R = reglas_espaciales(D, R)
    R = {k: v.fillna(False).astype(bool) for k, v in R.items()}
    print("Tarea C: métodos estadísticos y ML (C2)")
    M, Z, TEND, CPS = {}, {}, {}, []
    for k in ESTACIONES:
        fm, z, tend = anomalias_mad(D, k)
        mad = pd.DataFrame(fm).any(axis=1) if fm else pd.Series(False, index=D.idx)
        iflag, score = isolation_forest(D, k)
        reg, cp = regimenes(D, k)
        CPS.append(cp)
        ev = D.df[k]["presente"]
        M[k] = dict(reglas=reglas_sin_faltante(R[k]) & ev, mad=mad.astype(bool) & ev, iforest=iflag & ev, regimen=reg & ev)
        M[k]["estadistico"] = M[k]["mad"] | M[k]["iforest"] | M[k]["regimen"]
        M[k]["union"] = M[k]["reglas"] | M[k]["estadistico"]
        M[k]["etiqueta"] = D.etiqueta_debil(k)
        Z[k], TEND[k] = z, tend
    cps = pd.concat(CPS, ignore_index=True)
    cps.to_csv(os.path.join(out, "tabla_C2_cambios_de_regimen.csv"), index=False)
    print("Tarea C: comparación (C3)")
    met, acu, porc, desa = comparar(D, M)
    met.to_csv(os.path.join(out, "tabla_C3_metricas.csv"), index=False)
    acu.to_csv(os.path.join(out, "tabla_C3_acuerdo.csv"), index=False)
    porc.to_csv(os.path.join(out, "tabla_C3_deteccion_por_codigo.csv"), index=False)
    desa.sort_values(["estacion", "categoria", "minutos"], ascending=[True, True, False]).to_csv(
        os.path.join(out, "tabla_C3_desacuerdos.csv"), index=False)
    print("Tarea C: clasificación (C4)")
    cfg = hallazgos_configuracion(D)
    cfg.to_csv(os.path.join(out, "tabla_C4_hallazgos_configuracion.csv"), index=False)
    ev = clasificar_eventos(D, M, R, Z, TEND)
    ev.to_csv(os.path.join(out, "tabla_C4_eventos_clasificados.csv"), index=False)
    print("Tarea C: banderas y diccionario (C5)")
    B, textos, principal, bits = {}, {}, {}, {}
    for k in ESTACIONES:
        B[k], textos[k], principal[k], bits[k] = asignar_banderas(D, k, R[k], M[k], Z[k])
    BANDERAS.to_csv(os.path.join(out, "tabla_C5_esquema_banderas.csv"), index=False)
    pd.concat([inferir_codigos(D, k, R[k], M[k]) for k in ESTACIONES]).to_csv(
        os.path.join(out, "tabla_C5_diccionario_codigos_calidad.csv"), index=False)
    resumen = pd.DataFrame({k: principal[k].value_counts() for k in ESTACIONES}).fillna(0).astype(int)
    resumen.to_csv(os.path.join(out, "tabla_C5_resumen_banderas_principales.csv"))
    if exportar:
        for k in ESTACIONES:
            t = D.df[k][ESTACIONES[k]["vars"] + ["calidad"]].copy()
            t["banderas"], t["bandera_principal"], t["banderas_bits"], t["estado_dato"] = textos[k], principal[k], bits[k], "preliminar"
            t = t[D.df[k]["presente"]]
            t.index.name = "fecha_hora"
            t.to_csv(os.path.join(out, f"datos_banderas_{k}.csv.gz"), compression="gzip")
    fig_C1(D, R, out)
    fig_C2(D, M, Z, cps, out)
    fig_C3(met, acu, porc, out)
    fig_C4(ev, out)
    fig_C5(B, out)


# =====================================================================================
# 11. PROGRAMA PRINCIPAL
# =====================================================================================
def main():
    ap = argparse.ArgumentParser(description="EDA y detección de anomalías de las estaciones 35, 367, 803 y 6004")
    ap.add_argument("--datos", default=".", help="carpeta con los cuatro CSV")
    ap.add_argument("--salida", default="resultados", help="carpeta de salida")
    ap.add_argument("--solo-eda", action="store_true", help="ejecuta únicamente la Tarea A")
    ap.add_argument("--solo-anomalias", action="store_true", help="ejecuta únicamente la Tarea C")
    ap.add_argument("--exportar-datos", action="store_true", help="guarda los datos con banderas (csv.gz)")
    ap.add_argument("--lat", type=float, default=P["lat"], help="latitud de las estaciones (grados)")
    ap.add_argument("--lon", type=float, default=P["lon"], help="longitud de las estaciones (grados)")
    ap.add_argument("--desfase-6004", default="auto", help="horas a sumar a la marca de tiempo del 6004 (o 'auto')")
    a = ap.parse_args()
    P["lat"], P["lon"] = a.lat, a.lon
    os.makedirs(a.salida, exist_ok=True)
    D = Datos(a.datos, a.desfase_6004)
    print(f"Datos: {D.idx[0]} a {D.idx[-1]} | minutos esperados {len(D.idx)} | desfase 6004 {D.desfase_6004:+d} h")
    if not a.solo_anomalias:
        parte_A(D, a.salida)
    if not a.solo_eda:
        parte_C(D, a.salida, a.exportar_datos)
    print(f"Listo. Resultados en: {os.path.abspath(a.salida)}")


if __name__ == "__main__":
    main()
