import math
import tkinter as tk
from tkinter import ttk

# ======================== МАТРИЧНІ УТИЛІТИ 3×3 ========================
# Точка — рядок (x, y, 1), перетворення P' = P · M; у добутку A · B першим діє A

def mat_mul(A, B):
    C = [[0.0]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            C[i][j] = A[i][0]*B[0][j] + A[i][1]*B[1][j] + A[i][2]*B[2][j]
    return C

def T_shift(dx, dy):
    return [[1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [dx,  dy,  1.0]]

# Обертання навколо точки: T(-x0,-y0) · R(phi) · T(+x0,+y0)
def R_about(x0, y0, phi_deg):
    phi = math.radians(phi_deg)
    c, s = math.cos(phi), math.sin(phi)
    return mat_mul(
        mat_mul([[1, 0, 0], [0, 1, 0], [-x0, -y0, 1]],
                [[c, s, 0], [-s, c, 0], [0, 0, 1]]),
        [[1, 0, 0], [0, 1, 0], [x0, y0, 1]]
    )

# Симетрія відносно точки = обертання на 180° навколо неї
def S_point(px, py):
    return [[-1.0, 0.0, 0.0],
            [ 0.0,-1.0, 0.0],
            [ 2.0*px, 2.0*py, 1.0]]

# Масштабування відносно точки: T(-p) · diag(sx, sy, 1) · T(+p)
def S_about(px, py, sx, sy):
    S = [[sx, 0.0, 0.0],
         [0.0, sy, 0.0],
         [0.0, 0.0, 1.0]]
    return mat_mul(
        mat_mul([[1, 0, 0], [0, 1, 0], [-px, -py, 1]], S),
        [[1, 0, 0], [0, 1, 0], [px, py, 1]]
    )

def A_aff(a11, a12, a21, a22, tx, ty):
    return [[a11, a12, 0.0],
            [a21, a22, 0.0],
            [tx,  ty,  1.0]]

# Афінна матриця з образів осей і початку координат
def A_from_axes(Ox, Oy, VXx, VXy, VYx, VYy):
    return [[VXx, VYx, 0.0],
            [VXy, VYy, 0.0],
            [Ox,  Oy,  1.0]]

# P' = P · M з нормуванням на вагу W; W ≈ 0 — точка на лінії зникання
def apply33(M, x, y):
    W = x*M[0][2] + y*M[1][2] + 1*M[2][2]
    if abs(W) < 1e-12:
        return None
    X = (x*M[0][0] + y*M[1][0] + 1*M[2][0]) / W
    Y = (x*M[0][1] + y*M[1][1] + 1*M[2][1]) / W
    return (X, Y)

def fmt2(x):
    try:
        if isinstance(x, str):
            x = x.replace(",", ".")
        return f"{float(x):.2f}"
    except Exception:
        return "0.00"

# ======================== ГЕОМЕТРІЯ ДЕТАЛІ ========================
# Контур: 4 дуги з центрами в отворах + 2 трикутні шпилі угорі та внизу

def circle_points(cx, cy, r, steps=180):
    pts = []
    for i in range(steps + 1):
        ang = 2.0 * math.pi * i / steps
        pts.append((cx + r*math.cos(ang), cy + r*math.sin(ang)))
    return pts

class PartModel:
    def __init__(self):
        self.D_shell = 40.0
        self.D_hole  = 14.0
        self.tip_deg = 90.0
        self.hole_dx = 15.0
        self.hole_dy = 15.0
        self.x0_grid = 20.0
        self.y0_grid = 30.0
        self.nodes = {}
        self.segments = []

    def set_params(self, D_shell, D_hole, tip_deg, hole_dx, hole_dy):
        self.D_shell = float(D_shell)
        self.D_hole  = float(D_hole)
        self.tip_deg = float(tip_deg)
        self.hole_dx = float(hole_dx)
        self.hole_dy = float(hole_dy)

    def regenerate(self):
        self.nodes.clear()
        self.segments.clear()

        R  = max(1.0, self.D_shell / 2.0)
        r  = max(0.0, self.D_hole  / 2.0)
        hx = max(0.1, self.hole_dx)
        hy = max(0.1, self.hole_dy)

        # Висота шпиля за кутом при вершині
        deg = min(179.0, max(1.0, self.tip_deg))
        Ht  = hx / math.tan(math.radians(deg) / 2.0)

        # R > hy, інакше сусідні кола не перетнуться
        if R <= hy + 1e-3:
            R = hy + 1e-3

        sqy = math.sqrt(max(0.0, R*R - hy*hy))

        grid = 10.0
        cx_mod = round(self.x0_grid / grid) * grid + hx + sqy
        cy_mod = round(self.y0_grid / grid) * grid + hy + R + Ht

        UL = (cx_mod - hx, cy_mod + hy)
        UR = (cx_mod + hx, cy_mod + hy)
        BL = (cx_mod - hx, cy_mod - hy)
        BR = (cx_mod + hx, cy_mod - hy)

        UL_top = (cx_mod - hx, cy_mod + hy + R)
        UR_top = (cx_mod + hx, cy_mod + hy + R)
        BL_bot = (cx_mod - hx, cy_mod - hy - R)
        BR_bot = (cx_mod + hx, cy_mod - hy - R)

        T_apex = (cx_mod, cy_mod + hy + R + Ht)
        B_apex = (cx_mod, cy_mod - hy - R - Ht)

        L_cusp = (cx_mod - hx - sqy, cy_mod)
        R_cusp = (cx_mod + hx + sqy, cy_mod)

        self.nodes["TL"] = UL_top; self.nodes["TR"] = UR_top
        self.nodes["BL_pt"] = BL_bot; self.nodes["BR_pt"] = BR_bot
        self.nodes["TA"] = T_apex; self.nodes["BA"] = B_apex
        self.nodes["L"]  = L_cusp; self.nodes["R"]  = R_cusp
        self.nodes["UL"] = UL; self.nodes["UR"] = UR
        self.nodes["BL"] = BL; self.nodes["BR"] = BR
        self.center_model = (cx_mod, cy_mod)

        theta_y = math.atan2(hy, sqy)

        # Замкнений контур проти годинникової стрілки; [1:] прибирає дубль точки на стику
        outline = []
        outline.extend(self._arc_by_angles(UL, R, math.pi/2.0,        math.pi + theta_y,       steps=90))
        outline.extend(self._arc_by_angles(BL, R, math.pi - theta_y,  1.5*math.pi,             steps=90)[1:])
        outline.append(B_apex)
        outline.append(BR_bot)
        outline.extend(self._arc_by_angles(BR, R, -math.pi/2.0,       theta_y,                 steps=90)[1:])
        outline.extend(self._arc_by_angles(UR, R, -theta_y,           math.pi/2.0,             steps=90)[1:])
        outline.append(T_apex)
        outline.append(UL_top)
        self.segments.append(outline)

        if r > 0:
            for cx, cy in (UL, UR, BL, BR):
                self.segments.append(circle_points(cx, cy, r, steps=180))

    def _arc_by_angles(self, center, R, a1, a2, steps=180):
        cx, cy = center
        pts = []
        for i in range(steps + 1):
            t = i / steps
            a = a1 + (a2 - a1) * t
            pts.append((cx + R*math.cos(a), cy + R*math.sin(a)))
        return pts


# ==================== ДОДАТок ========================

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Valeriia Rozumovych-Lab1 variant 16")
        self.geometry("1280x840")

        self.model = PartModel()
        self.model.regenerate()

        # Підсумкова матриця: M_obj · M_scale · (M_world або M_proj)
        self.M_obj   = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        self.M_scale = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        self.M_world = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

        self._density   = tk.IntVar(value=35)   # px/cm; 1 см = 10 одиниць
        self.show_grid  = tk.BooleanVar(value=True)

        self.flag_eu = tk.BooleanVar(value=False)
        self.flag_af = tk.BooleanVar(value=False)
        self.flag_pr = tk.BooleanVar(value=False)

        c = self.model.center_model
        self.center_xy = [c[0], c[1]]
        self.rotate_xy = [c[0] + 60.0, c[1] + 40.0]
        self.eu_fig_x = tk.DoubleVar(value=self.center_xy[0])
        self.eu_fig_y = tk.DoubleVar(value=self.center_xy[1])
        self.eu_rot_x = tk.DoubleVar(value=self.rotate_xy[0])
        self.eu_rot_y = tk.DoubleVar(value=self.rotate_xy[1])
        self.eu_angle = tk.DoubleVar(value=0.0)

        self.sym_pt = [c[0], c[1]]
        self.eu_sym_x = tk.DoubleVar(value=c[0])
        self.eu_sym_y = tk.DoubleVar(value=c[1])

        self.drag_target = None

        self.aff_xX = tk.StringVar(value="1.00")
        self.aff_xY = tk.StringVar(value="0.00")
        self.aff_yX = tk.StringVar(value="0.00")
        self.aff_yY = tk.StringVar(value="1.00")
        self.aff_x0 = tk.StringVar(value="0.00")
        self.aff_y0 = tk.StringVar(value="0.00")

        self.scale_x = tk.DoubleVar(value=1.0)
        self.scale_y = tk.DoubleVar(value=1.0)

        # Ручки афінного базису довжиною 3 клітинки
        self._aff_cell = 10.0
        self._aff_handle_cells = 3.0
        self._nx = 1; self._ny = 1
        self.af_O = [0.0, 0.0]
        self.af_X = [30.0, 0.0]
        self.af_Y = [0.0, 30.0]

        self._pr_cell = 10.0
        self.pr_xx = tk.DoubleVar(value=18.0)
        self.pr_xy = tk.DoubleVar(value=0.0)
        self.pr_yx = tk.DoubleVar(value=0.0)
        self.pr_yy = tk.DoubleVar(value=18.0)
        self.pr_x0 = tk.DoubleVar(value=0.0)
        self.pr_y0 = tk.DoubleVar(value=0.0)
        self.pr_wX = tk.DoubleVar(value=1.0)
        self.pr_wY = tk.DoubleVar(value=1.0)
        self.pr_w0 = tk.DoubleVar(value=350.0)

        self.pr_O = [self.pr_x0.get()/max(1e-9, self.pr_w0.get()),
                     self.pr_y0.get()/max(1e-9, self.pr_w0.get())]
        self.pr_X = [self.pr_xx.get()*self._pr_cell/max(1e-9, self.pr_wX.get()), 0.0]
        self.pr_Y = [0.0, self.pr_yy.get()*self._pr_cell/max(1e-9, self.pr_wY.get())]

        self._axis_gap_px   = 18
        self._mask_overlap  = 4

        self._build_ui()
        self._bind_canvas()
        self._sync_aff_ui_from_Mworld()
        self._sync_aff_handles_from_Mworld()
        self._redraw()

    # ---------- Права панель із параметрами ----------
    def _build_ui(self):
        paned = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True)

        self.canvas = tk.Canvas(paned, bg="white")
        paned.add(self.canvas, weight=4)

        right = ttk.Frame(paned, width=430); right.pack_propagate(False)
        paned.add(right, weight=1)

        sc = tk.Canvas(right, highlightthickness=0)
        vs = ttk.Scrollbar(right, orient="vertical", command=sc.yview)
        wrap = ttk.Frame(sc)
        wrap.bind("<Configure>", lambda e: sc.configure(scrollregion=sc.bbox("all")))
        sc.create_window((0, 0), window=wrap, anchor="nw", width=410)
        sc.configure(yscrollcommand=vs.set); sc.pack(side="left", fill="both", expand=True)
        vs.pack(side="right", fill="y")

        dims = ttk.LabelFrame(wrap, text="Figure size & Display"); dims.pack(fill=tk.X, padx=10, pady=8)
        self._add_dim(dims, "Shell arc Ø:",  "_od", self.model.D_shell)
        self._add_dim(dims, "Hole Ø:",       "_id", self.model.D_hole)
        self._add_dim(dims, "Tip angle:", "_ta", self.model.tip_deg)
        self._add_dim(dims, "Hole dx:",      "_hx", self.model.hole_dx)
        self._add_dim(dims, "Hole dy:",      "_hy", self.model.hole_dy)

        scale_frame = ttk.LabelFrame(dims, text="Scale (about figure center)")
        scale_frame.pack(fill=tk.X, padx=6, pady=6)
        row_sc = ttk.Frame(scale_frame); row_sc.pack(fill=tk.X, padx=4, pady=4)
        ttk.Label(row_sc, text="Sx").pack(side=tk.LEFT)
        ttk.Entry(row_sc, textvariable=self.scale_x, width=7).pack(side=tk.LEFT, padx=(6, 14))
        ttk.Label(row_sc, text="Sy").pack(side=tk.LEFT)
        ttk.Entry(row_sc, textvariable=self.scale_y, width=7).pack(side=tk.LEFT, padx=6)
        ttk.Button(row_sc, text="Apply scale", command=self._apply_scale).pack(side=tk.RIGHT)
        ttk.Button(row_sc, text="Reset scale", command=self._reset_scale).pack(side=tk.RIGHT, padx=(0, 6))

        row = ttk.Frame(dims); row.pack(fill=tk.X, padx=6, pady=4)
        ttk.Label(row, text="Density (px/cm):").pack(side=tk.LEFT)
        ttk.Entry(row, textvariable=self._density, width=6).pack(side=tk.LEFT, padx=6)
        self._density.trace_add("write", lambda *_: self._redraw())
        ttk.Button(row, text="Set default figure", command=self._set_default_figure).pack(side=tk.RIGHT)

        flags = ttk.LabelFrame(wrap, text="Transformation"); flags.pack(fill=tk.X, padx=10, pady=8)
        ttk.Checkbutton(flags, text="Euclidean",  variable=self.flag_eu, command=self._toggle_panels).pack(anchor="w")
        ttk.Checkbutton(flags, text="Affine",     variable=self.flag_af, command=self._toggle_panels).pack(anchor="w")
        ttk.Checkbutton(flags, text="Projective", variable=self.flag_pr, command=self._toggle_panels).pack(anchor="w")

        self.panel_eu = ttk.LabelFrame(wrap, text="Euclidean parameters")
        r1 = ttk.Frame(self.panel_eu); r1.pack(fill=tk.X, padx=6, pady=3)
        self._xy_row(r1, "Figure position:", self.eu_fig_x, self.eu_fig_y, self._eu_apply_pos)
        r2 = ttk.Frame(self.panel_eu); r2.pack(fill=tk.X, padx=6, pady=3)
        self._xy_row(r2, "Rotate point position:", self.eu_rot_x, self.eu_rot_y, self._eu_apply_rotpos)
        ttk.Button(self.panel_eu, text="Reset Euclidean", command=self._reset_eu).pack(fill=tk.X, padx=6, pady=4)
        r3 = ttk.Frame(self.panel_eu); r3.pack(fill=tk.X, padx=6, pady=3)
        ttk.Label(r3, text="Angle (deg):").pack(side=tk.LEFT)
        ttk.Entry(r3, textvariable=self.eu_angle, width=6).pack(side=tk.LEFT, padx=6)
        ttk.Button(r3, text="Apply rotation", command=self._eu_apply_rotation).pack(side=tk.RIGHT)
        r4 = ttk.Frame(self.panel_eu); r4.pack(fill=tk.X, padx=6, pady=6)
        ttk.Label(r4, text="Symmetry point:").pack(side=tk.LEFT)
        ttk.Label(r4, text="X").pack(side=tk.LEFT, padx=(8, 2))
        ttk.Entry(r4, textvariable=self.eu_sym_x, width=7).pack(side=tk.LEFT)
        ttk.Label(r4, text="Y").pack(side=tk.LEFT, padx=(8, 2))
        ttk.Entry(r4, textvariable=self.eu_sym_y, width=7).pack(side=tk.LEFT)
        ttk.Button(r4, text="Apply symmetry", command=self._eu_apply_symmetry).pack(side=tk.RIGHT)

        self.panel_af = ttk.LabelFrame(wrap, text="Affine parameters")
        grid = ttk.Frame(self.panel_af); grid.pack(fill=tk.X, padx=6, pady=4)
        self._cell(grid, "xX:", self.aff_xX, 0, 0); self._cell(grid, "yX:", self.aff_yX, 0, 2); self._const_cell(grid, "0", 0, 4)
        self._cell(grid, "xY:", self.aff_xY, 1, 0); self._cell(grid, "yY:", self.aff_yY, 1, 2); self._const_cell(grid, "0", 1, 4)
        self._cell(grid, "x0:", self.aff_x0, 2, 0); self._cell(grid, "y0:", self.aff_y0, 2, 2); self._const_cell(grid, "1", 2, 4)
        row_af = ttk.Frame(self.panel_af); row_af.pack(fill=tk.X, padx=6, pady=4)
        ttk.Button(row_af, text="Apply affine", command=self._af_apply_set).pack(side=tk.LEFT, expand=True, fill=tk.X)
        ttk.Button(row_af, text="Reset affine", command=self._af_reset).pack(side=tk.RIGHT, expand=True, fill=tk.X)

        self.panel_pr = ttk.LabelFrame(wrap, text="Projective parameters")
        gridp = ttk.Frame(self.panel_pr); gridp.pack(fill=tk.X, padx=6, pady=4)
        self._cell(gridp, "xx:", self.pr_xx, 0, 0); self._cell(gridp, "yx:", self.pr_yx, 0, 2); self._cell(gridp, "wX:", self.pr_wX, 0, 4)
        self._cell(gridp, "xy:", self.pr_xy, 1, 0); self._cell(gridp, "yy:", self.pr_yy, 1, 2); self._cell(gridp, "wY:", self.pr_wY, 1, 4)
        self._cell(gridp, "x0:", self.pr_x0, 2, 0); self._cell(gridp, "y0:", self.pr_y0, 2, 2); self._cell(gridp, "w0:", self.pr_w0, 2, 4)
        row_pr = ttk.Frame(self.panel_pr); row_pr.pack(fill=tk.X, padx=6, pady=4)
        ttk.Button(row_pr, text="Apply projective", command=self._pr_apply_from_ui).pack(side=tk.LEFT, expand=True, fill=tk.X)
        ttk.Button(row_pr, text="Reset projective", command=self._pr_reset).pack(side=tk.RIGHT, expand=True, fill=tk.X)

        self._toggle_panels()

    def _add_dim(self, parent, label, name, init):
        row = ttk.Frame(parent); row.pack(fill=tk.X, padx=6, pady=3)
        ttk.Label(row, text=label).pack(side=tk.LEFT)
        var = tk.DoubleVar(value=float(init))
        ttk.Entry(row, textvariable=var, width=10).pack(side=tk.RIGHT)
        setattr(self, name, var)
        var.trace_add("write", lambda *_: self._on_dims_change())

    def _xy_row(self, parent, label, xv, yv, cb):
        ttk.Label(parent, text=label).pack(side=tk.LEFT)
        ttk.Label(parent, text="X").pack(side=tk.LEFT, padx=(8, 2))
        ttk.Entry(parent, textvariable=xv, width=7).pack(side=tk.LEFT)
        ttk.Label(parent, text="Y").pack(side=tk.LEFT, padx=(8, 2))
        ttk.Entry(parent, textvariable=yv, width=7).pack(side=tk.LEFT)
        ttk.Button(parent, text="Apply", command=cb).pack(side=tk.RIGHT)

    def _cell(self, parent, text, var, r, c):
        f = ttk.Frame(parent); f.grid(row=r, column=c, padx=4, pady=3, sticky="w")
        ttk.Label(f, text=text).pack(side=tk.LEFT)
        e = ttk.Entry(f, textvariable=var, width=8)
        e.pack(side=tk.LEFT)
        def _normalize(_evt=None, _v=var):
            _v.set(fmt2(_v.get()))
        e.bind("<FocusOut>", _normalize)

    def _const_cell(self, parent, value, r, c):
        f = ttk.Frame(parent); f.grid(row=r, column=c, padx=4, pady=3, sticky="w")
        e = ttk.Entry(f, width=6); e.insert(0, value); e.state(["disabled"]); e.pack(side=tk.LEFT)

    def _toggle_panels(self):
        for panel, flag in ((self.panel_eu, self.flag_eu),
                            (self.panel_af, self.flag_af),
                            (self.panel_pr, self.flag_pr)):
            if flag.get(): panel.pack(fill=tk.X, padx=10, pady=4)
            else:          panel.pack_forget()
        self._redraw()

    # ---------- Події полотна ----------
    def _bind_canvas(self):
        self.canvas.bind("<Configure>", lambda e: self._redraw())
        self.canvas.bind("<MouseWheel>", self._on_wheel)
        self.canvas.bind("<Button-4>",  self._on_wheel)
        self.canvas.bind("<Button-5>",  self._on_wheel)
        self.canvas.bind("<ButtonPress-1>", self._on_press)
        self.canvas.bind("<B1-Motion>",     self._on_drag)
        self.canvas.bind("<ButtonRelease-1>", lambda e: setattr(self, "drag_target", None))

    def _on_dims_change(self):
        try:
            self.model.set_params(
                D_shell=float(self._od.get()),
                D_hole =float(self._id.get()),
                tip_deg=float(self._ta.get()),
                hole_dx=float(self._hx.get()),
                hole_dy=float(self._hy.get()),
            )
        except Exception:
            return
        self.model.regenerate()
        self._redraw()

    def _set_default_figure(self):
        self._od.set(40.0)
        self._id.set(14.0)
        self._ta.set(90.0)
        self._hx.set(15.0)
        self._hy.set(15.0)

    # Масштабування відносно поточного центра фігури
    def _apply_scale(self):
        try:
            sx = float(self.scale_x.get())
            sy = float(self.scale_y.get())
        except Exception:
            return
        c_world = apply33(self.M_obj, *self.model.center_model)
        if c_world:
            px, py = c_world
            self.M_scale = mat_mul(self.M_scale, S_about(px, py, sx, sy))
            self._redraw()

    def _reset_scale(self):
        self.M_scale = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        self.scale_x.set(1.0)
        self.scale_y.set(1.0)
        self._redraw()

    # ---------- Евклідові перетворення ----------
    def _eu_apply_pos(self):
        target_x = self.eu_fig_x.get(); target_y = self.eu_fig_y.get()
        cur = apply33(self.M_obj, *self.model.center_model)
        if cur:
            dx, dy = target_x - cur[0], target_y - cur[1]
            self.M_obj = mat_mul(self.M_obj, T_shift(dx, dy))
        self._redraw()

    def _eu_apply_rotpos(self):
        self.rotate_xy[:] = [self.eu_rot_x.get(), self.eu_rot_y.get()]
        self._redraw()

    def _eu_apply_rotation(self):
        self.M_obj = mat_mul(self.M_obj, R_about(self.rotate_xy[0], self.rotate_xy[1], self.eu_angle.get()))
        self._redraw()

    def _eu_apply_symmetry(self):
        self.sym_pt[:] = [float(self.eu_sym_x.get()), float(self.eu_sym_y.get())]
        px, py = self.sym_pt
        self.M_obj = mat_mul(self.M_obj, S_point(px, py))
        self._redraw()

    def _reset_eu(self):
        c = self.model.center_model
        self.center_xy[:] = [c[0], c[1]]
        self.rotate_xy[:] = [c[0] + 60.0, c[1] + 40.0]
        self.eu_fig_x.set(c[0]); self.eu_fig_y.set(c[1])
        self.eu_rot_x.set(self.rotate_xy[0]); self.eu_rot_y.set(self.rotate_xy[1])
        self.eu_angle.set(0.0)
        self.sym_pt[:] = [c[0], c[1]]
        self.eu_sym_x.set(c[0]); self.eu_sym_y.set(c[1])
        self.M_obj = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        self._redraw()

    # ---------- Афінні перетворення: UI ↔ матриця ----------
    # У полях UI значення задані в клітинках (1 клітинка = 10 одиниць)
    def _ui_to_world_matrix(self):
        u = self._aff_cell
        a11 = float(self.aff_xX.get()) / u
        a12 = float(self.aff_xY.get()) / u
        a21 = float(self.aff_yX.get()) / u
        a22 = float(self.aff_yY.get()) / u
        tx  = float(self.aff_x0.get()) * u
        ty  = float(self.aff_y0.get()) * u
        return A_aff(a11, a12, a21, a22, tx, ty)

    def _world_matrix_to_ui(self):
        u = self._aff_cell
        self.aff_xX.set(fmt2(self.M_world[0][0] * u))
        self.aff_xY.set(fmt2(self.M_world[0][1] * u))
        self.aff_yX.set(fmt2(self.M_world[1][0] * u))
        self.aff_yY.set(fmt2(self.M_world[1][1] * u))
        self.aff_x0.set(fmt2(self.M_world[2][0] / u))
        self.aff_y0.set(fmt2(self.M_world[2][1] / u))

    def _af_apply_set(self):
        self.M_world = self._ui_to_world_matrix()
        self._sync_aff_handles_from_Mworld()
        self._redraw()

    def _af_reset(self):
        self.M_world = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        self._sync_aff_ui_from_Mworld()
        self._sync_aff_handles_from_Mworld()
        self._redraw()

    # ---------- Світ → екран ----------
    def _grid_metrics(self):
        w = self.canvas.winfo_width(); h = self.canvas.winfo_height()
        pad = getattr(self, "grid_padding_px", 48)
        left, right  = pad, w - pad
        top, bottom  = pad, h - pad
        s = max(1, int(self._density.get())) / 10.0
        step_world = 10.0
        step_exact = step_world * s
        return left, right, top, bottom, step_exact, s, step_world

    # Вісь Y у Canvas спрямована вниз, тому інвертуємо
    def _world_to_screen(self, X, Y):
        left, right, top, bottom, step_exact, s, step_world = self._grid_metrics()
        return (round(left + X * s), round(bottom - Y * s))

    # ======================== ПРОЄКТИВНИЙ РЕЖИМ ========================
    # Ваги (wX, wY, w0) у третьому стовпці роблять перетворення проєктивним
    def _M_proj(self):
        k = self._pr_cell
        return [[self.pr_xx.get()*k*self.pr_wX.get(), self.pr_xy.get()*k*self.pr_wX.get(), self.pr_wX.get()],
                [self.pr_yx.get()*k*self.pr_wY.get(), self.pr_yy.get()*k*self.pr_wY.get(), self.pr_wY.get()],
                [self.pr_x0.get()*k*self.pr_w0.get(), self.pr_y0.get()*k*self.pr_w0.get(), self.pr_w0.get()]]

    def _pr_apply_from_ui(self):
        k = self._pr_cell
        w0 = max(1e-9, self.pr_w0.get())
        self.pr_O[:] = [self.pr_x0.get()/w0, self.pr_y0.get()/w0]
        self.pr_X[:] = [self.pr_xx.get()*k/max(1e-9, self.pr_wX.get()), 0.0]
        self.pr_Y[:] = [0.0, self.pr_yy.get()*k/max(1e-9, self.pr_wY.get())]
        self._redraw()

    def _pr_reset(self):
        self.pr_xx.set(18.0); self.pr_yy.set(18.0)
        self.pr_xy.set(0.0);  self.pr_yx.set(0.0)
        self.pr_x0.set(0.0);  self.pr_y0.set(0.0)
        self.pr_wX.set(1.0);  self.pr_wY.set(1.0); self.pr_w0.set(350.0)
        self._pr_apply_from_ui()

    # Вісь обрізається на 98 % шляху до точки зникання
    def _axis_u_end(self):
        xx = self.pr_xx.get()*self._pr_cell; x0 = self.pr_x0.get()
        wX = self.pr_wX.get(); w0 = self.pr_w0.get()
        if abs(wX) < 1e-9:
            return 400.0
        x_v = xx / wX
        x_t = x_v * 0.98
        denom = (xx - x_t*wX)
        if abs(denom) < 1e-9: return 400.0
        u = (x_t*w0 - x0) / denom
        return max(40.0, min(5000.0, u))

    def _axis_v_end(self):
        yy = self.pr_yy.get()*self._pr_cell; y0 = self.pr_y0.get()
        wY = self.pr_wY.get(); w0 = self.pr_w0.get()
        if abs(wY) < 1e-9:
            return 400.0
        y_v = yy / wY
        y_t = y_v * 0.98
        denom = (yy - y_t*wY)
        if abs(denom) < 1e-9: return 400.0
        v = (y_t*w0 - y0) / denom
        return max(40.0, min(5000.0, v))

    # ---------- Сітка, осі, ручки ----------
    # Сітка малюється як образ системи координат під поточним перетворенням
    def _draw_grid(self):
        if not self.show_grid.get(): return

        if not self.flag_pr.get():
            left, right, top, bottom, step_exact, s, step_world = self._grid_metrics()
            nx = max(1, math.floor((right - left) / step_exact))
            ny = max(1, math.floor((bottom - top) / step_exact))
            self._nx, self._ny = nx, ny

            def F(u, v): return apply33(self.M_world, u, v)

            for i in range(0, nx + 1):
                u = i * step_world
                x0, y0 = F(u, 0.0); x1, y1 = F(u, ny * step_world)
                sx0, sy0 = self._world_to_screen(x0, y0); sx1, sy1 = self._world_to_screen(x1, y1)
                self.canvas.create_line(sx0, sy0, sx1, sy1, fill="#ededed")

            for j in range(0, ny + 1):
                v = j * step_world
                x0, y0 = F(0.0, v); x1, y1 = F(nx * step_world, v)
                sx0, sy0 = self._world_to_screen(x0, y0); sx1, sy1 = self._world_to_screen(x1, y1)
                self.canvas.create_line(sx0, sy0, sx1, sy1, fill="#ededed")

            X0 = F(0.0, 0.0); X1 = F(nx * step_world, 0.0)
            sx0, sy0 = self._world_to_screen(*X0); sx1, sy1 = self._world_to_screen(*X1)
            self.canvas.create_line(sx0, sy0, sx1, sy1, fill="red", width=2, arrow=tk.LAST)
            self.canvas.create_text(sx1-10, sy1+14, text="X", fill="red")

            Y0 = F(0.0, 0.0); Y1 = F(0.0, ny * step_world)
            s0x, s0y = self._world_to_screen(*Y0); s1x, s1y = self._world_to_screen(*Y1)
            self.canvas.create_line(s0x, s0y, s1x, s1y, fill="green", width=2, arrow=tk.LAST)
            self.canvas.create_text(s1x-12, s1y-10, text="Y", fill="green")

            if self.flag_af.get():
                self._sync_aff_handles_from_Mworld()
                for (pt, color) in ((self.af_O, "#222"), (self.af_X, "#b91c1c"), (self.af_Y, "#15803d")):
                    sx, sy = self._world_to_screen(*pt)
                    self.canvas.create_oval(sx-5, sy-5, sx+5, sy+5, fill=color, outline="")

            if self.flag_eu.get():
                c_world = apply33(self.M_obj, *self.model.center_model)
                if c_world:
                    cx, cy = self._world_to_screen(*c_world)
                    self.canvas.create_oval(cx-5, cy-5, cx+5, cy+5, fill="#18c5de", outline="")
                rx, ry = self._world_to_screen(*self.rotate_xy)
                self.canvas.create_oval(rx-5, ry-5, rx+5, ry+5, fill="#ffd84d", outline="")
                sx, sy = self._world_to_screen(*self.sym_pt)
                self.canvas.create_oval(sx-5, sy-5, sx+5, sy+5, fill="#e11d48", outline="")
            return

        left, right, top, bottom, step_exact, s, step_world = self._grid_metrics()
        Umax = self._axis_u_end()
        Vmax = self._axis_v_end()
        M = self._M_proj()
        F = lambda u, v: apply33(M, u, v)

        u_lines = max(1, int(Umax/step_world))
        v_lines = max(1, int(Vmax/step_world))
        for i in range(u_lines+1):
            u = i*step_world; p0, p1 = F(u, 0.0), F(u, Vmax)
            if p0 and p1:
                x0, y0 = self._world_to_screen(*p0); x1, y1 = self._world_to_screen(*p1)
                self.canvas.create_line(x0, y0, x1, y1, fill="#ededed")
        for j in range(v_lines+1):
            v = j*step_world; p0, p1 = F(0.0, v), F(Umax, v)
            if p0 and p1:
                x0, y0 = self._world_to_screen(*p0); x1, y1 = self._world_to_screen(*p1)
                self.canvas.create_line(x0, y0, x1, y1, fill="#ededed")

        O  = F(0.0, 0.0); X1 = F(Umax, 0.0); Y1 = F(0.0, Vmax)

        def arrow_seg(p0, p1, gap):
            sx0, sy0 = self._world_to_screen(*p0); sx1, sy1 = self._world_to_screen(*p1)
            dx, dy = sx1-sx0, sy1-sy0; L = max(1.0, math.hypot(dx, dy))
            ex, ey = sx1 - dx/L*gap, sy1 - dy/L*gap
            return (sx0, sy0, ex, ey), (dx, dy)

        if O and X1:
            (sx0, sy0, ex, ey), (dx, dy) = arrow_seg(O, X1, self._axis_gap_px)
            self.canvas.create_line(sx0, sy0, ex, ey, fill="red", width=2, arrow=tk.LAST)
            self.canvas.create_text(ex-10, ey+14, text="X", fill="red")
            # Біла маска: затирає сітку за наконечником стрілки
            nxp, nyp = -dy, dx; K = 6000
            self.canvas.create_polygon(ex+self._mask_overlap, ey,
                                       ex+dx, ey+dy,
                                       ex+dx+nxp*K, ey+dy+nyp*K,
                                       ex+nxp*K, ey+nyp*K, fill="white", outline="")

        if O and Y1:
            (sx0, sy0, ex, ey), (dx, dy) = arrow_seg(O, Y1, self._axis_gap_px)
            self.canvas.create_line(sx0, sy0, ex, ey, fill="green", width=2, arrow=tk.LAST)
            self.canvas.create_text(ex-12, ey-10, text="Y", fill="green")
            nxp, nyp = dy, -dx; K = 6000
            self.canvas.create_polygon(ex, ey+self._mask_overlap,
                                       ex+dx, ey+dy,
                                       ex+dx+nxp*K, ey+dy+nyp*K,
                                       ex+nxp*K, ey+nyp*K, fill="white", outline="")

        self.pr_O[:] = [self.pr_x0.get()/max(1e-9, self.pr_w0.get()),
                        self.pr_y0.get()/max(1e-9, self.pr_w0.get())]
        self.pr_X[:] = [self.pr_xx.get()*self._pr_cell/max(1e-9, self.pr_wX.get()), 0.0]
        self.pr_Y[:] = [0.0, self.pr_yy.get()*self._pr_cell/max(1e-9, self.pr_wY.get())]
        for pt in (self.pr_O, self.pr_Y, self.pr_X):
            sx, sy = self._world_to_screen(*pt)
            self.canvas.create_oval(sx-6, sy-6, sx+6, sy+6, fill="#111", outline="")

    # Евклід  масштаб  афін проєктив
    def _redraw(self):
        self.canvas.delete("all")
        self._draw_grid()
        M = mat_mul(self.M_obj, self.M_scale)
        if self.flag_pr.get():
            M = mat_mul(M, self._M_proj())
        else:
            M = mat_mul(M, self.M_world)
        for poly in self.model.segments:
            pts = []
            for (x, y) in poly:
                P = apply33(M, x, y)
                if P is None: continue
                sx, sy = self._world_to_screen(*P)
                pts.extend([sx, sy])
            if len(pts) >= 4:
                self.canvas.create_line(*pts, fill="#203a8a", width=2,
                                        capstyle=tk.ROUND, joinstyle=tk.ROUND)

    def _on_wheel(self, e):
        factor = 1.1 if (getattr(e, "delta", 0) > 0 or getattr(e, "num", 0) == 4) else 1/1.1
        cur = max(1.0, float(self._density.get()))
        new_density = max(1.0, min(1000.0, cur * factor))
        self._density.set(int(round(new_density)))
        self._redraw()

    def _from_screen_to_world(self, sx, sy):
        left, right, top, bottom, step_exact, s, step_world = self._grid_metrics()
        return ((sx - left) / s, (bottom - sy) / s)

    # ---------- Перетягування ручок ----------
    def _on_press(self, e):
        x, y = e.x, e.y
        def hit_world(pt, r=12):
            sx, sy = self._world_to_screen(*pt)
            return (sx-x)**2 + (sy-y)**2 <= r*r
        self.drag_target = None

        if self.flag_eu.get():
            c_world = apply33(self.M_obj, *self.model.center_model)
            if c_world and hit_world(c_world): self.drag_target = "center"; return
            if hit_world(self.rotate_xy):      self.drag_target = "rot";    return
            if hit_world(self.sym_pt):         self.drag_target = "sym";    return

        if self.flag_af.get():
            if hit_world(self.af_O): self.drag_target = "aff_O"; return
            if hit_world(self.af_X): self.drag_target = "aff_X"; return
            if hit_world(self.af_Y): self.drag_target = "aff_Y"; return

        if self.flag_pr.get():
            if hit_world(self.pr_O): self.drag_target = "pr_O"; return
            if hit_world(self.pr_X): self.drag_target = "pr_X"; return
            if hit_world(self.pr_Y): self.drag_target = "pr_Y"; return

    def _on_drag(self, e):
        if not self.drag_target: return
        Xw, Yw = self._from_screen_to_world(e.x, e.y)

        if self.drag_target == "center":
            cur = apply33(self.M_obj, *self.model.center_model)
            if cur:
                dx, dy = Xw - cur[0], Yw - cur[1]
                self.M_obj = mat_mul(self.M_obj, T_shift(dx, dy))
                new_center = apply33(self.M_obj, *self.model.center_model)
                if new_center:
                    self.eu_fig_x.set(round(new_center[0], 3))
                    self.eu_fig_y.set(round(new_center[1], 3))

        elif self.drag_target == "rot":
            self.rotate_xy[:] = [Xw, Yw]
            self.eu_rot_x.set(Xw); self.eu_rot_y.set(Yw)

        elif self.drag_target == "sym":
            self.sym_pt[:] = [Xw, Yw]
            self.eu_sym_x.set(round(Xw, 3)); self.eu_sym_y.set(round(Yw, 3))

        elif self.drag_target in ("aff_O", "aff_X", "aff_Y"):
            VXx, VXy = self.M_world[0][0], self.M_world[1][0]
            VYx, VYy = self.M_world[0][1], self.M_world[1][1]
            Ox,  Oy  = self.M_world[2][0], self.M_world[2][1]

            Lx = self._aff_handle_cells * self._aff_cell
            Ly = self._aff_handle_cells * self._aff_cell

            if self.drag_target == "aff_O":
                Ox, Oy = Xw, Yw
            elif self.drag_target == "aff_X":
                VXx = (Xw - Ox) / Lx
                VXy = (Yw - Oy) / Lx
            else:
                VYx = (Xw - Ox) / Ly
                VYy = (Yw - Oy) / Ly

            self.M_world = A_from_axes(Ox, Oy, VXx, VXy, VYx, VYy)
            self._world_matrix_to_ui()

        elif self.drag_target == "pr_O":
            self.pr_O[:] = [Xw, Yw]
            self.pr_x0.set(Xw * self.pr_w0.get())
            self.pr_y0.set(Yw * self.pr_w0.get())

        elif self.drag_target == "pr_X":
            # Положення ручки задає вагу wX
            Xw = max(1e-6, Xw)
            self.pr_X[:] = [Xw, 0.0]
            self.pr_wX.set((self.pr_xx.get()*self._pr_cell) / Xw)

        elif self.drag_target == "pr_Y":
            Yw = max(1e-6, Yw)
            self.pr_Y[:] = [0.0, Yw]
            self.pr_wY.set((self.pr_yy.get()*self._pr_cell) / Yw)

        self._sync_aff_handles_from_Mworld()
        self._redraw()

    def _sync_aff_ui_from_Mworld(self):
        self._world_matrix_to_ui()

    def _sync_aff_handles_from_Mworld(self):
        Lx = self._aff_handle_cells * self._aff_cell
        Ly = self._aff_handle_cells * self._aff_cell
        self.af_O[:] = apply33(self.M_world, 0.0, 0.0)
        self.af_X[:] = apply33(self.M_world, Lx, 0.0)
        self.af_Y[:] = apply33(self.M_world, 0.0, Ly)


if __name__ == "__main__":
    App().mainloop()
