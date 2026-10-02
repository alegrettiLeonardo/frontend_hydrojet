import os
import subprocess
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from datetime import datetime


# ----------------------------------------------------------------------
# CONFIGURAÇÃO DO BINÁRIO DO HYDROJET
# ----------------------------------------------------------------------
HYDROJET_BIN_DIR = "/home/leonardo/Documentos/hydrojet/bin"
HYDROJET_EXE_NAME = "hydrojet"   # ajuste aqui se o binário tiver outro nome
HYDROJET_DEFAULT_NAME = "DEFAULT.txt"


# ----------------------------------------------------------------------
# Definição dos parâmetros (nome, rótulo, grupo, tipo, default, choices)
# ----------------------------------------------------------------------

def now_datetime_str():
    return datetime.now().strftime("%d/%m/%Y. %H:%M:%S")


PARAMS = [
    # ----- Cabeçalho -----
    dict(name="TITLE", label="Título do caso", group="Geral",
         ptype="text", default="1 pad (5 pocket) Default data water HJB"),
    dict(name="DATETIME", label="Data/hora", group="Geral",
         ptype="text", default=now_datetime_str()),

    # ----- Linha 3 -----
    dict(name="NPADS", label="NPADS (nº de pads)", group="Geral",
         ptype="int", default=2),
    dict(name="NPOCKETS_PAD", label="NPOCKETS_PAD (pockets por pad)", group="Geral",
         ptype="int", default=0),
    dict(name="NXT", label="NXT (pontos circunferenciais)", group="Geral",
         ptype="int", default=25),
    dict(name="NPC", label="NPC (pontos recess circunf.)", group="Geral",
         ptype="int", default=0),
    dict(name="NLA_NYI", label="NLA/NYI (pontos axiais land)", group="Geral",
         ptype="int", default=8),
    dict(name="NPA", label="NPA (pontos axiais recess)", group="Geral",
         ptype="int", default=0),

    dict(name="IFULL", label="IFULL (1=360°,0=pads)", group="Geral",
         ptype="choice", default="0 - Pads",
         choices=["0 - Pads", "1 - 360°"]),
    dict(name="ITAPER", label="ITAPER (folga cônica)", group="Geral",
         ptype="choice", default="2 - Tapered padrão",
         choices=["0 - Sem taper", "1 - Outro", "2 - Tapered padrão"]),
    dict(name="IROW_HST", label="IROW_HST (fileira hidrostat.)", group="Geral",
         ptype="choice", default="2 - Single row",
         choices=["1 - Outra", "2 - Single row"]),
    dict(name="NLOBES", label="NLOBES (lóbulos selo)", group="Geral",
         ptype="int", default=1),
    dict(name="ITYPE_PAD", label="ITYPE_PAD (1=tilting,0=fixed)", group="Geral",
         ptype="choice", default="0 - Fixed pad",
         choices=["0 - Fixed pad", "1 - Tilting pad"]),

    # ----- Linha 4 -----
    dict(name="IFILM_INERTIA", label="IFILM_INERTIA (inércia no filme)", group="Geral",
         ptype="choice", default="1 - Com inércia",
         choices=["0 - Sem inércia", "1 - Com inércia"]),
    dict(name="IINLET_INERTIA", label="IINLET_INERTIA (inércia entrada/saída)", group="Geral",
         ptype="choice", default="0 - Sem inércia",
         choices=["0 - Sem inércia", "1 - Com inércia"]),
    dict(name="MAXIT_LAND", label="MAXIT_LAND (máx. iterações land)", group="Numérico",
         ptype="int", default=999),
    dict(name="MAXIT_RECESS", label="MAXIT_RECESS (máx. iterações recess)", group="Numérico",
         ptype="int", default=10),

    dict(name="FLUID_TYPE", label="FLUID_TYPE (tipo de fluido)", group="Fluido",
         ptype="choice", default="7 - Óleo",
         choices=[
             "1 - H2", "2 - N2", "3 - O2", "4 - CH4",
             "5 - Ar", "6 - Água", "7 - Óleo", "8 - Outro"
         ]),
    dict(name="IWORN", label="IWORN (mancal desgastado)", group="Geral",
         ptype="choice", default="0 - Sem desgaste",
         choices=["0 - Sem desgaste", "1 - Desgastado"]),
    dict(name="ISTEPPED_CLEAR", label="ISTEPPED_CLEAR (folga escalonada)", group="Geral",
         ptype="choice", default="0 - Sem degrau",
         choices=["0 - Sem degrau", "1 - Com degrau"]),
    dict(name="BEARING_TYPE", label="BEARING_TYPE (1=hidrostat.,2=selo,3=hidrodin.)", group="Geral",
         ptype="choice", default="3 - Hidrodinâmico",
         choices=["1 - Hidrostático", "2 - Selo anular", "3 - Hidrodinâmico"]),
    dict(name="THERMAL_OPTION", label="THERMAL_OPTION (modelo térmico)", group="Fluido",
         ptype="choice", default="1 - Isotérmico",
         choices=[
             "-1 - Especial", "0 - Adiabático", "1 - Isotérmico",
             "2 - Misto 2", "3 - Misto 3", "4 - Radial 1", "5 - Radial 2"
         ]),
    dict(name="LIFT_FLAG1", label="LIFT_FLAG1", group="Numérico",
         ptype="choice", default="1 - Ativo",
         choices=["0 - Inativo", "1 - Ativo"]),
    dict(name="LIFT_FLAG2", label="LIFT_FLAG2", group="Numérico",
         ptype="choice", default="0 - Inativo",
         choices=["0 - Inativo", "1 - Ativo"]),

    # ----- Linha 5 -----
    dict(name="IPRUNI", label="IPRUNI (pressão uniforme direita)", group="Numérico",
         ptype="choice", default="1 - Uniforme",
         choices=["0 - Fourier", "1 - Uniforme"]),
    dict(name="IPLUNI", label="IPLUNI (pressão uniforme esquerda)", group="Numérico",
         ptype="choice", default="1 - Uniforme",
         choices=["0 - Fourier", "1 - Uniforme"]),
    dict(name="NPRCS", label="NPRCS (# Fourier direita)", group="Numérico",
         ptype="int", default=0),
    dict(name="NPLCS", label="NPLCS (# Fourier esquerda)", group="Numérico",
         ptype="int", default=0),

    # ----- Linha 6 – Geometria -----
    dict(name="CINLET", label="CINLET (folga entrada)", group="Geometria",
         ptype="float", default=0.00010425),
    dict(name="CEXIT", label="CEXIT (folga saída)", group="Geometria",
         ptype="float", default=0.00010425),
    dict(name="D_JOURNAL", label="D_JOURNAL (diâmetro rotor)", group="Geometria",
         ptype="float", default=0.125),
    dict(name="L_TOTAL", label="L_TOTAL (comprimento axial)", group="Geometria",
         ptype="float", default=0.085),
    dict(name="AR", label="AR (comprimento recess)", group="Geometria",
         ptype="float", default=0.0),
    dict(name="HREC", label="HREC (profundidade recess)", group="Geometria",
         ptype="float", default=0.0),
    dict(name="VSUP_ORIF", label="VSUP_ORIF (volume orifício)", group="Geometria",
         ptype="float", default=0.0),
    dict(name="LR", label="LR (comprimento à direita)", group="Geometria",
         ptype="float", default=0.0425),

    # ----- Pads -----
    dict(name="PAD1_NREC", label="Pad1 NREC", group="Geometria",
         ptype="int", default=0),
    dict(name="PAD1_LEAD", label="Pad1 leading (°)", group="Geometria",
         ptype="float", default=24.0),
    dict(name="PAD1_ARC", label="Pad1 arco (°)", group="Geometria",
         ptype="float", default=132.0),
    dict(name="PAD1_ROT", label="Pad1 rot (rad)", group="Geometria",
         ptype="float", default=0.0),

    dict(name="PAD2_NREC", label="Pad2 NREC", group="Geometria",
         ptype="int", default=0),
    dict(name="PAD2_LEAD", label="Pad2 leading (°)", group="Geometria",
         ptype="float", default=204.0),
    dict(name="PAD2_ARC", label="Pad2 arco (°)", group="Geometria",
         ptype="float", default=132.0),
    dict(name="PAD2_ROT", label="Pad2 rot (rad)", group="Geometria",
         ptype="float", default=0.0),

    # ----- Alinhamento (simplificado) -----
    dict(name="EX", label="EX (excentricidade X)", group="Geometria",
         ptype="float", default=1e-20),
    dict(name="EY", label="EY (excentricidade Y)", group="Geometria",
         ptype="float", default=0.00010425),
    dict(name="GEO_AUX1", label="GEO_AUX1", group="Geometria",
         ptype="float", default=0.085),
    dict(name="GEO_AUX2", label="GEO_AUX2", group="Geometria",
         ptype="float", default=0.00010425),

    dict(name="ALIGN1", label="ALIGN1", group="Geometria",
         ptype="float", default=0.0),
    dict(name="ALIGN2", label="ALIGN2", group="Geometria",
         ptype="float", default=0.0),
    dict(name="ALIGN3", label="ALIGN3", group="Geometria",
         ptype="float", default=0.0),
    dict(name="ALIGN4", label="ALIGN4", group="Geometria",
         ptype="float", default=0.0),
    dict(name="ALIGN5", label="ALIGN5", group="Geometria",
         ptype="float", default=0.0),
    dict(name="ALIGN6", label="ALIGN6", group="Geometria",
         ptype="float", default=1e-20),
    dict(name="ALIGN7", label="ALIGN7", group="Geometria",
         ptype="float", default=1e-20),

    # ----- Propriedades fluido -----
    dict(name="MU_SUP", label="MU_SUP (μ supr.)", group="Fluido",
         ptype="float", default=0.01694),
    dict(name="RHO_SUP", label="RHO_SUP (ρ supr.)", group="Fluido",
         ptype="float", default=876.0),
    dict(name="MU_DISC", label="MU_DISC (μ desc.)", group="Fluido",
         ptype="float", default=0.01694),
    dict(name="RHO_DISC", label="RHO_DISC (ρ desc.)", group="Fluido",
         ptype="float", default=876.0),
    dict(name="COMPRESS", label="COMPRESS", group="Fluido",
         ptype="float", default=0.0000000006),

    dict(name="CP", label="CP", group="Fluido",
         ptype="float", default=1959.0),
    dict(name="K_FLUID", label="K_FLUID", group="Fluido",
         ptype="float", default=0.131),
    dict(name="BETA_T", label="BETA_T", group="Fluido",
         ptype="float", default=0.0008),
    dict(name="VTC", label="VTC", group="Fluido",
         ptype="float", default=0.0102),

    # ----- Operação -----
    dict(name="RPM", label="RPM", group="Operação",
         ptype="float", default=99000.0),

    dict(name="P1", label="P1 (P_supply)", group="Operação",
         ptype="float", default=100000.0),
    dict(name="RHO1", label="RHO1", group="Operação",
         ptype="float", default=876.0),
    dict(name="MU1", label="MU1", group="Operação",
         ptype="float", default=0.01694),
    dict(name="KLOSS1", label="KLOSS1", group="Operação",
         ptype="float", default=0.0),

    dict(name="P2", label="P2", group="Operação",
         ptype="float", default=100000.0),
    dict(name="RHO2", label="RHO2", group="Operação",
         ptype="float", default=876.0),
    dict(name="MU2", label="MU2", group="Operação",
         ptype="float", default=0.01694),
    dict(name="KLOSS2", label="KLOSS2", group="Operação",
         ptype="float", default=0.0),

    dict(name="P3_SUP", label="P3_SUP", group="Operação",
         ptype="float", default=100000.0),
    dict(name="P3_DISCH", label="P3_DISCH", group="Operação",
         ptype="float", default=100000.0),
    dict(name="MU3_SUP", label="MU3_SUP", group="Operação",
         ptype="float", default=0.01694),
    dict(name="MU3_DISCH", label="MU3_DISCH", group="Operação",
         ptype="float", default=0.01694),
    dict(name="RHO3_SUP", label="RHO3_SUP", group="Operação",
         ptype="float", default=876.0),
    dict(name="RHO3_DISCH", label="RHO3_DISCH", group="Operação",
         ptype="float", default=876.0),

    # ----- Orifício / swirl -----
    dict(name="CD_ORIF", label="CD_ORIF", group="Operação",
         ptype="float", default=1.0),
    dict(name="D_ORIF", label="D_ORIF", group="Operação",
         ptype="float", default=0.0),
    dict(name="K_INLET_SEAL", label="K_INLET_SEAL", group="Operação",
         ptype="float", default=0.2),
    dict(name="SWIRL_INLET", label="SWIRL_INLET", group="Operação",
         ptype="float", default=0.5),
    dict(name="P_REC_RATIO", label="P_REC_RATIO", group="Operação",
         ptype="float", default=0.5),

    # ----- Perdas pocket / bordas -----
    dict(name="K_POCKET_UP", label="K_POCKET_UP", group="Operação",
         ptype="float", default=0.0),
    dict(name="K_POCKET_DOWN", label="K_POCKET_DOWN", group="Operação",
         ptype="float", default=0.0),
    dict(name="K_POCKET_AXIAL", label="K_POCKET_AXIAL", group="Operação",
         ptype="float", default=0.0),
    dict(name="K_EXTRA_SEAL", label="K_EXTRA_SEAL", group="Operação",
         ptype="float", default=0.0),
    dict(name="K_LEADING_EDGE", label="K_LEADING_EDGE", group="Operação",
         ptype="float", default=0.0),

    # ----- Turbulência / rugosidade -----
    dict(name="T_SUPPLY", label="T_SUPPLY (K)", group="Fluido",
         ptype="float", default=333.9),
    dict(name="ROUGH_ROTOR", label="ROUGH_ROTOR", group="Fluido",
         ptype="float", default=1e-20),
    dict(name="ROUGH_STATOR", label="ROUGH_STATOR", group="Fluido",
         ptype="float", default=0.0),
    dict(name="A_TURB", label="A_TURB", group="Fluido",
         ptype="float", default=0.0),
    dict(name="B_TURB", label="B_TURB", group="Fluido",
         ptype="float", default=0.001375),
    dict(name="TURB_SCALE", label="TURB_SCALE", group="Fluido",
         ptype="float", default=500000.0),
    dict(name="EXPO_TURB", label="EXPO_TURB", group="Fluido",
         ptype="float", default=0.3333),

    # ----- Relaxação / erros -----
    dict(name="RELAX_MOM", label="RELAX_MOM", group="Numérico",
         ptype="float", default=0.8),
    dict(name="RELAX_P", label="RELAX_P", group="Numérico",
         ptype="float", default=0.6),
    dict(name="RELAX_T", label="RELAX_T", group="Numérico",
         ptype="float", default=0.9),
    dict(name="ERR_MASS_LAND", label="ERR_MASS_LAND", group="Numérico",
         ptype="float", default=0.01),
    dict(name="ERR_P_RECESS", label="ERR_P_RECESS", group="Numérico",
         ptype="float", default=0.01),
    dict(name="ERR_EXTRA", label="ERR_EXTRA", group="Numérico",
         ptype="float", default=0.0001),

    # ----- Preload / mistura / complacência -----
    dict(name="PRELOAD", label="PRELOAD", group="Geometria",
         ptype="float", default=0.0),
    dict(name="OFFSET_PIVOT", label="OFFSET_PIVOT (0-1)", group="Geometria",
         ptype="float", default=0.5),
    dict(name="THERM_MIX", label="THERM_MIX", group="Fluido",
         ptype="float", default=0.5),
    dict(name="COMP_P_FACTOR", label="COMP_P_FACTOR", group="Numérico",
         ptype="float", default=0.0),
    dict(name="COMP_LOSS_FACTOR", label="COMP_LOSS_FACTOR", group="Numérico",
         ptype="float", default=0.0),
    dict(name="COMP_RELAX", label="COMP_RELAX", group="Numérico",
         ptype="float", default=1.0),

    # ----- T rotor / pad / honeycomb -----
    dict(name="T_ROTOR", label="T_ROTOR", group="Fluido",
         ptype="float", default=363.3),
    dict(name="T_PAD_INNER", label="T_PAD_INNER", group="Fluido",
         ptype="float", default=353.3),
    dict(name="HONEYCOMB_DEPTH", label="HONEYCOMB_DEPTH", group="Geometria",
         ptype="float", default=0.0),

    # ----- Casos dinâmicos (N_CASES, FREQ_OPTION) -----
    dict(name="N_CASES", label="N_CASES", group="Casos",
         ptype="int", default=9),
    dict(name="FREQ_OPTION", label="FREQ_OPTION (1=síncrona,2=não)", group="Casos",
         ptype="choice", default="2 - Não síncrona",
         choices=["1 - Síncrona", "2 - Não síncrona"]),
]


# Acrescenta CASE1..CASE9
for i, rpm_default in enumerate([400, 850, 1300, 1750, 2200, 2650, 3100, 3550, 4000], start=1):
    PARAMS.extend([
        dict(name=f"CASE{i}_SPEED", label=f"Caso {i} - RPM", group="Casos",
             ptype="float", default=float(rpm_default)),
        dict(name=f"CASE{i}_FLAG1", label=f"Caso {i} - FLAG1", group="Casos",
             ptype="choice", default="1 - Ativo",
             choices=["0 - Inativo", "1 - Ativo"]),
        dict(name=f"CASE{i}_FLAG2", label=f"Caso {i} - FLAG2", group="Casos",
             ptype="choice", default="1 - Ativo",
             choices=["0 - Inativo", "1 - Ativo"]),
        dict(name=f"CASE{i}_FLAG3", label=f"Caso {i} - FLAG3", group="Casos",
             ptype="choice", default="0 - Inativo",
             choices=["0 - Inativo", "1 - Ativo"]),
        dict(name=f"CASE{i}_LIMIT", label=f"Caso {i} - LIMIT", group="Casos",
             ptype="float", default=10000.0),
    ])


# ----------------------------------------------------------------------
# Classe da GUI
# ----------------------------------------------------------------------

class HydrojetGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Hydrojet DEFAULT.txt - GUI")
        self.geometry("1100x700")

        self.param_specs = {p["name"]: p for p in PARAMS}
        self.param_vars = {}

        # Vars da aba Execução
        self.bin_dir_var = tk.StringVar(value=HYDROJET_BIN_DIR)
        self.exe_name_var = tk.StringVar(value=HYDROJET_EXE_NAME)
        self.exec_log = None  # será criado na aba Execução

        self._build_ui()

    # --------------- Construção da interface -----------------

    def _build_ui(self):
        notebook = ttk.Notebook(self)
        notebook.pack(fill=tk.BOTH, expand=True)

        # Criar abas
        tabs = {}
        for tab_name in [
            "Geral", "Geometria", "Fluido", "Operação",
            "Numérico", "Casos", "Execução", "Prévia DEFAULT.txt"
        ]:
            frame = ttk.Frame(notebook)
            notebook.add(frame, text=tab_name)
            tabs[tab_name] = frame

        # Frames por grupo
        group_to_tab = {
            "Geral": "Geral",
            "Geometria": "Geometria",
            "Fluido": "Fluido",
            "Operação": "Operação",
            "Numérico": "Numérico",
            "Casos": "Casos",
        }

        # Criar widgets para cada parâmetro
        grid_pos = {key: 0 for key in group_to_tab}  # row por grupo

        for spec in PARAMS:
            name = spec["name"]
            group = spec["group"]
            tab_key = group_to_tab.get(group, "Geral")
            frame = tabs[tab_key]

            row = grid_pos[group]
            grid_pos[group] += 1

            label = ttk.Label(frame, text=spec["label"])
            label.grid(row=row, column=0, sticky="w", padx=4, pady=2)

            var = tk.StringVar()
            # default como string
            var.set(str(spec["default"]))
            self.param_vars[name] = var

            if spec["ptype"] == "choice" and "choices" in spec:
                combo = ttk.Combobox(frame, textvariable=var, state="readonly")
                combo["values"] = spec["choices"]
                combo.grid(row=row, column=1, sticky="we", padx=4, pady=2)
            else:
                entry = ttk.Entry(frame, textvariable=var)
                entry.grid(row=row, column=1, sticky="we", padx=4, pady=2)

            # expandir coluna 1
            frame.grid_columnconfigure(1, weight=1)

        # --------------- Aba de Execução -----------------
        exec_frame = tabs["Execução"]

        row = 0
        ttk.Label(exec_frame, text="Pasta do hydrojet:").grid(
            row=row, column=0, sticky="w", padx=4, pady=2
        )
        entry_bin = ttk.Entry(exec_frame, textvariable=self.bin_dir_var)
        entry_bin.grid(row=row, column=1, sticky="we", padx=4, pady=2)
        btn_browse_bin = ttk.Button(exec_frame, text="Procurar...", command=self.browse_bin_dir)
        btn_browse_bin.grid(row=row, column=2, sticky="w", padx=4, pady=2)

        row += 1
        ttk.Label(exec_frame, text="Nome do executável:").grid(
            row=row, column=0, sticky="w", padx=4, pady=2
        )
        entry_exe = ttk.Entry(exec_frame, textvariable=self.exe_name_var)
        entry_exe.grid(row=row, column=1, sticky="we", padx=4, pady=2)

        exec_frame.grid_columnconfigure(1, weight=1)

        row += 1
        btns_exec = ttk.Frame(exec_frame)
        btns_exec.grid(row=row, column=0, columnspan=3, sticky="we", padx=4, pady=4)

        btn_save_bin = ttk.Button(btns_exec, text="Salvar DEFAULT na pasta", command=self.save_default_to_bin)
        btn_save_bin.pack(side=tk.LEFT, padx=4)

        btn_run = ttk.Button(btns_exec, text="Rodar Hydrojet", command=self.run_hydrojet_only)
        btn_run.pack(side=tk.LEFT, padx=4)

        btn_save_run = ttk.Button(btns_exec, text="Salvar & rodar Hydrojet", command=self.save_and_run_hydrojet)
        btn_save_run.pack(side=tk.LEFT, padx=4)

        row += 1
        ttk.Label(exec_frame, text="Log de execução:").grid(
            row=row, column=0, columnspan=3, sticky="w", padx=4, pady=(8, 2)
        )

        row += 1
        self.exec_log = tk.Text(exec_frame, wrap="none", height=15)
        self.exec_log.grid(row=row, column=0, columnspan=3, sticky="nsew", padx=4, pady=4)

        yscroll_exec = ttk.Scrollbar(exec_frame, orient=tk.VERTICAL, command=self.exec_log.yview)
        yscroll_exec.grid(row=row, column=3, sticky="ns")
        self.exec_log.configure(yscrollcommand=yscroll_exec.set)

        xscroll_exec = ttk.Scrollbar(exec_frame, orient=tk.HORIZONTAL, command=self.exec_log.xview)
        xscroll_exec.grid(row=row + 1, column=0, columnspan=3, sticky="we")
        self.exec_log.configure(xscrollcommand=xscroll_exec.set)

        exec_frame.grid_rowconfigure(row, weight=1)
        exec_frame.grid_columnconfigure(1, weight=1)

        # --------------- Aba de prévia -----------------
        preview_frame = tabs["Prévia DEFAULT.txt"]
        self.preview_text = tk.Text(preview_frame, wrap="none")
        self.preview_text.pack(fill=tk.BOTH, expand=True)

        # Barra de rolagem horizontal/vertical na prévia
        yscroll = ttk.Scrollbar(preview_frame, orient=tk.VERTICAL, command=self.preview_text.yview)
        yscroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.preview_text.configure(yscrollcommand=yscroll.set)

        xscroll = ttk.Scrollbar(preview_frame, orient=tk.HORIZONTAL, command=self.preview_text.xview)
        xscroll.pack(side=tk.BOTTOM, fill=tk.X)
        self.preview_text.configure(xscrollcommand=xscroll.set)

        # Botões embaixo (globais)
        btn_frame = ttk.Frame(self)
        btn_frame.pack(fill=tk.X, pady=4)

        btn_preview = ttk.Button(btn_frame, text="Atualizar prévia", command=self.update_preview)
        btn_preview.pack(side=tk.LEFT, padx=4)

        btn_save = ttk.Button(btn_frame, text="Salvar DEFAULT.txt", command=self.save_default_file)
        btn_save.pack(side=tk.LEFT, padx=4)

        # usa a mesma função da aba Execução
        btn_run_all = ttk.Button(btn_frame, text="Salvar & rodar Hydrojet", command=self.save_and_run_hydrojet)
        btn_run_all.pack(side=tk.LEFT, padx=4)

        btn_exit = ttk.Button(btn_frame, text="Sair", command=self.destroy)
        btn_exit.pack(side=tk.RIGHT, padx=4)

        # Prévia inicial
        self.update_preview()

    # --------------- Utilitários -----------------

    def _get_value_str(self, name):
        """Retorna o valor em string já preparado para o DEFAULT.txt."""
        spec = self.param_specs[name]
        raw = self.param_vars[name].get().strip()

        # Para choices, pegar só o código numérico (antes do espaço)
        if spec["ptype"] == "choice":
            if " " in raw:
                raw = raw.split(" ", 1)[0]

        # Números: trocar vírgula por ponto
        if spec["ptype"] in ("int", "float"):
            raw = raw.replace(",", ".")

        # Não fazer mais nada, Fortran lê em formato livre
        return raw

    def browse_bin_dir(self):
        new_dir = filedialog.askdirectory(
            initialdir=self.bin_dir_var.get() or HYDROJET_BIN_DIR,
            title="Selecionar pasta do hydrojet"
        )
        if new_dir:
            self.bin_dir_var.set(new_dir)

    def _append_exec_log(self, text):
        if self.exec_log is None:
            return
        self.exec_log.insert(tk.END, text)
        self.exec_log.see(tk.END)

    # --------------- Geração do DEFAULT.txt -----------------

    def build_default_lines(self):
        g = self._get_value_str  # atalho

        lines = []

        # Linha 1 - título
        lines.append(g("TITLE"))
        # Linha 2 - data/hora
        lines.append(g("DATETIME"))

        # Linha 3
        lines.append(" ".join([
            g("NPADS"), g("NPOCKETS_PAD"), g("NXT"), g("NPC"),
            g("NLA_NYI"), g("NPA"),
            g("IFULL"), g("ITAPER"), g("IROW_HST"), g("NLOBES"), g("ITYPE_PAD")
        ]))

        # Linha 4
        lines.append(" ".join([
            g("IFILM_INERTIA"), g("IINLET_INERTIA"), g("MAXIT_LAND"),
            g("MAXIT_RECESS"), g("FLUID_TYPE"), g("IWORN"),
            g("ISTEPPED_CLEAR"), g("BEARING_TYPE"), g("THERMAL_OPTION"),
            g("LIFT_FLAG1"), g("LIFT_FLAG2")
        ]))

        # Linha 5
        lines.append(" ".join([
            g("IPRUNI"), g("IPLUNI"), g("NPRCS"), g("NPLCS")
        ]))

        # Linha 6
        lines.append(" ".join([
            g("CINLET"), g("CEXIT"), g("D_JOURNAL"), g("L_TOTAL"),
            g("AR"), g("HREC"), g("VSUP_ORIF"), g("LR")
        ]))

        # Linhas 7–8 pads
        lines.append(" ".join([
            g("PAD1_NREC"), g("PAD1_LEAD"), g("PAD1_ARC"), g("PAD1_ROT")
        ]))
        lines.append(" ".join([
            g("PAD2_NREC"), g("PAD2_LEAD"), g("PAD2_ARC"), g("PAD2_ROT")
        ]))

        # Linhas 9–11 alinhamento/aux
        lines.append(" ".join([g("EX"), g("EY")]))
        lines.append(" ".join([g("GEO_AUX1"), g("GEO_AUX2")]))
        lines.append(" ".join([
            g("ALIGN1"), g("ALIGN2"), g("ALIGN3"),
            g("ALIGN4"), g("ALIGN5"), g("ALIGN6"), g("ALIGN7")
        ]))

        # Linhas 12–13 propriedades fluido
        lines.append(" ".join([
            g("MU_SUP"), g("RHO_SUP"), g("MU_DISC"),
            g("RHO_DISC"), g("COMPRESS")
        ]))
        lines.append(" ".join([
            g("CP"), g("K_FLUID"), g("BETA_T"), g("VTC")
        ]))

        # Linha 14 – RPM
        lines.append(g("RPM"))

        # Linhas 15–17 – pressões extras
        lines.append(" ".join([g("P1"), g("RHO1"), g("MU1"), g("KLOSS1")]))
        lines.append(" ".join([g("P2"), g("RHO2"), g("MU2"), g("KLOSS2")]))
        lines.append(" ".join([
            g("P3_SUP"), g("P3_DISCH"), g("MU3_SUP"),
            g("MU3_DISCH"), g("RHO3_SUP"), g("RHO3_DISCH")
        ]))

        # Linha 18 – orifício/swirl
        lines.append(" ".join([
            g("CD_ORIF"), g("D_ORIF"), g("K_INLET_SEAL"),
            g("SWIRL_INLET"), g("P_REC_RATIO")
        ]))

        # Linha 19 – perdas pocket/bordas
        lines.append(" ".join([
            g("K_POCKET_UP"), g("K_POCKET_DOWN"),
            g("K_POCKET_AXIAL"), g("K_EXTRA_SEAL"), g("K_LEADING_EDGE")
        ]))

        # Linha 20 – T, rugosidade, turbulência
        lines.append(" ".join([
            g("T_SUPPLY"), g("ROUGH_ROTOR"), g("ROUGH_STATOR"),
            g("A_TURB"), g("B_TURB"), g("TURB_SCALE"), g("EXPO_TURB")
        ]))

        # Linha 21 – relaxação/erros
        lines.append(" ".join([
            g("RELAX_MOM"), g("RELAX_P"), g("RELAX_T"),
            g("ERR_MASS_LAND"), g("ERR_P_RECESS"), g("ERR_EXTRA")
        ]))

        # Linha 22 – preload/offset
        lines.append(" ".join([g("PRELOAD"), g("OFFSET_PIVOT")]))

        # Linha 23 – mistura térmica
        lines.append(g("THERM_MIX"))

        # Linha 24 – complacência
        lines.append(" ".join([
            g("COMP_P_FACTOR"), g("COMP_LOSS_FACTOR"), g("COMP_RELAX")
        ]))

        # Linha 25 – temperaturas
        lines.append(" ".join([g("T_ROTOR"), g("T_PAD_INNER")]))

        # Linha 26 – honeycomb
        lines.append(g("HONEYCOMB_DEPTH"))

        # Linha 27 – N_CASES, FREQ_OPTION
        lines.append(" ".join([
            g("N_CASES"), g("FREQ_OPTION")
        ]))

        # Linhas 28–(28+N_CASES-1) – casos dinâmicos
        n_cases_str = g("N_CASES")
        try:
            n_cases = int(float(n_cases_str))
        except ValueError:
            n_cases = 0

        for i in range(1, n_cases + 1):
            lines.append(" ".join([
                g(f"CASE{i}_SPEED"),
                g(f"CASE{i}_FLAG1"),
                g(f"CASE{i}_FLAG2"),
                g(f"CASE{i}_FLAG3"),
                g(f"CASE{i}_LIMIT"),
            ]))

        return lines

    # --------------- Callbacks -----------------

    def update_preview(self):
        lines = self.build_default_lines()
        text = "\n".join(lines)
        self.preview_text.delete("1.0", tk.END)
        self.preview_text.insert("1.0", text)

    def save_default_file(self):
        """Salvar DEFAULT.txt em local escolhido pelo usuário."""
        lines = self.build_default_lines()
        text = "\n".join(lines)

        path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            title="Salvar DEFAULT.txt"
        )
        if not path:
            return

        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write(text)
            messagebox.showinfo("OK", f"Arquivo salvo em:\n{path}")
        except Exception as e:
            messagebox.showerror("Erro ao salvar", str(e))

    def save_default_to_bin(self):
        """Salvar DEFAULT.txt diretamente na pasta do hydrojet (aba Execução)."""
        lines = self.build_default_lines()
        text = "\n".join(lines)

        bin_dir = self.bin_dir_var.get().strip() or HYDROJET_BIN_DIR

        if not os.path.isdir(bin_dir):
            messagebox.showerror("Erro", f"Pasta do hydrojet não encontrada:\n{bin_dir}")
            return

        default_path = os.path.join(bin_dir, HYDROJET_DEFAULT_NAME)

        try:
            with open(default_path, "w", encoding="utf-8") as f:
                f.write(text)
            messagebox.showinfo("OK", f"DEFAULT.txt salvo em:\n{default_path}")
        except Exception as e:
            messagebox.showerror("Erro ao salvar DEFAULT.txt", str(e))

    def save_and_run_hydrojet(self):
        """Salvar DEFAULT.txt na pasta configurada e rodar o hydrojet."""
        self.save_default_to_bin()
        # mesmo se der erro, tentaremos rodar; poderia checar com try/except,
        # mas o caminho será validado novamente em _run_hydrojet
        self.run_hydrojet_only()

    def run_hydrojet_only(self):
        """Rodar o hydrojet usando DEFAULT.txt já existente na pasta."""
        bin_dir = self.bin_dir_var.get().strip() or HYDROJET_BIN_DIR
        exe_name = self.exe_name_var.get().strip() or HYDROJET_EXE_NAME
        default_path = os.path.join(bin_dir, HYDROJET_DEFAULT_NAME)

        if not os.path.isdir(bin_dir):
            messagebox.showerror("Erro", f"Pasta do hydrojet não encontrada:\n{bin_dir}")
            return

        if not os.path.isfile(default_path):
            messagebox.showerror(
                "Erro",
                f"DEFAULT.txt não encontrado em:\n{default_path}\n"
                f"Salve o arquivo antes de executar."
            )
            return

        self._run_hydrojet(bin_dir, exe_name)

    def _run_hydrojet(self, bin_dir, exe_name):
        exe_path = os.path.join(bin_dir, exe_name)

        if not os.path.isfile(exe_path):
            messagebox.showerror(
                "Erro",
                f"Binário do hydrojet não encontrado:\n{exe_path}\n"
                f"Verifique o nome em 'Nome do executável'."
            )
            return

        # Limpa log antes de uma nova execução
        if self.exec_log is not None:
            self.exec_log.delete("1.0", tk.END)

        try:
            result = subprocess.run(
                [exe_path],
                cwd=bin_dir,
                capture_output=True,
                text=True
            )
        except Exception as e:
            messagebox.showerror("Erro ao executar hydrojet", str(e))
            self._append_exec_log(f"Erro ao executar: {e}\n")
            return

        # Escreve log na aba Execução
        self._append_exec_log(f"Comando: {exe_path}\n")
        self._append_exec_log(f"CWD: {bin_dir}\n")
        self._append_exec_log(f"Retorno: {result.returncode}\n\n")

        if result.stdout:
            self._append_exec_log("----- STDOUT -----\n")
            self._append_exec_log(result.stdout + "\n")
        if result.stderr:
            self._append_exec_log("----- STDERR -----\n")
            self._append_exec_log(result.stderr + "\n")

        if result.returncode == 0:
            messagebox.showinfo(
                "Hydrojet",
                "Execução concluída com sucesso.\n"
                "Veja o log na aba 'Execução'."
            )
        else:
            messagebox.showwarning(
                "Hydrojet",
                f"Hydrojet terminou com código {result.returncode}.\n"
                f"Veja o log na aba 'Execução' para detalhes."
            )


if __name__ == "__main__":
    app = HydrojetGUI()
    app.mainloop()