import math
import tkinter as tk
from tkinter import ttk

TWO_PI = 2.0 * math.pi

DEFAULT_A = 40.0
DEFAULT_H = 60.0
DEFAULT_T_MIN = -4.0 * math.pi
DEFAULT_T_MAX = 4.0 * math.pi
DEFAULT_DT = 0.02
DEFAULT_T0 = 2.0
MAX_POINTS = 8000
COLORS = ["black", "red", "orange", "green", "blue", "indigo", "violet"]


# ---------------- чисельні методи ----------------

# інтеграл методом сімпсона: кінці з вагою 1, внутрішні точки по черзі з вагами 4 і 2
def simpson(f, a, b, n):
    if n % 2:
        n += 1
    step = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4 if i % 2 else 2) * f(a + i * step)
    return s * step / 3.0


# повний еліптичний інтеграл другого роду E від параметра m через середнє арифметико-геометричне
def elliptic_e(m):
    if m >= 1.0:
        return 1.0
    a, b = 1.0, math.sqrt(1.0 - m)
    power = 0.5
    s = power * m
    for _ in range(60):
        if abs(a - b) <= 1e-15 * a:
            break
        c = (a - b) / 2.0
        a, b = (a + b) / 2.0, math.sqrt(a * b)
        power *= 2.0
        s += power * c * c
    return math.pi / (2.0 * a) * (1.0 - s)


# ---------------- модель трохоїди ----------------

# колесо радіуса a котиться по осі X, точка жорстко прив'язана на відстані h від центра
class Trochoid:
    def __init__(self, a=DEFAULT_A, h=DEFAULT_H):
        self.a = float(a)
        self.h = float(h)

    # параметричне рівняння: x = a·t − h·sin t, y = a − h·cos t
    def point(self, t):
        return (self.a * t - self.h * math.sin(t),
                self.a - self.h * math.cos(t))

    # перша похідна, вектор швидкості: x' = a − h·cos t, y' = h·sin t
    def d1(self, t):
        return self.a - self.h * math.cos(t), self.h * math.sin(t)

    # друга похідна: x'' = h·sin t, y'' = h·cos t
    def d2(self, t):
        return self.h * math.sin(t), self.h * math.cos(t)

    # модуль швидкості
    def speed(self, t):
        dx, dy = self.d1(t)
        return math.hypot(dx, dy)

    # вістря- швидкість нульова, буває лише у циклоїди при t = 2πk
    def is_cusp(self, t):
        return self.speed(t) < 1e-6 * max(1.0, self.a + self.h)

    # одиничний вектор дотичної, у вістрі граничний напрям вектор другої похідної
    def unit_tangent(self, t):
        dx, dy = self.d2(t) if self.is_cusp(t) else self.d1(t)
        norm = math.hypot(dx, dy)
        return dx / norm, dy / norm

    # векторний добуток x'y'' − y'x'', після підстановки h·a·cos t − h²
    def cross(self, t):
        return self.h * (self.a * math.cos(t) - self.h)

    # ρ = |r'|³ / |x'y'' − y'x''|, у вістрі нуль, при нульовій кривині нескінченність
    def radius_of_curvature(self, t):
        if self.is_cusp(t):
            return 0.0
        denominator = abs(self.cross(t))
        if denominator < 1e-12:
            return math.inf
        return self.speed(t) ** 3 / denominator

    # вкорочена трохоїда: точка всередині колеса
    def is_curtate(self):
        return 0.0 < self.h < self.a - 1e-9

    # циклоїда: точка на ободі
    def is_cycloid(self):
        return abs(self.h - self.a) <= 1e-9

    # подовжена трохоїда: точка за ободом, є петлі
    def is_prolate(self):
        return self.h > self.a + 1e-9

    # тип кривої для поля «тип»
    def kind(self):
        if self.h == 0:
            return "пряма y = a, h = 0"
        if self.is_cycloid():
            return "циклоїда, h = a"
        if self.is_curtate():
            return "вкорочена трохоїда, h < a"
        return "подовжена трохоїда, h > a"

    # точки перегину: cos t = h/a, тобто t = ±arccos h/a + 2πk, лише при 0 < h < a
    def inflections(self, t1, t2):
        if not self.is_curtate():
            return []
        base = math.acos(self.h / self.a)
        result = []
        for t0 in (base, -base):
            k = math.ceil((t1 - t0) / TWO_PI)
            while t0 + k * TWO_PI <= t2:
                result.append(t0 + k * TWO_PI)
                k += 1
        result.sort()
        return result

    # довжина дуги як інтеграл модуля швидкості методом сімпсона
    def arc_length(self, t1, t2):
        n = min(20000, max(1000, int((t2 - t1) * 200)))
        return simpson(self.speed, t1, t2, n)

    # точна довжина за період: L = 4·S·E від m, де S = a + h, m = 4ah / S²
    def period_length_exact(self):
        total = self.a + self.h
        return 4.0 * total * elliptic_e(4.0 * self.a * self.h / (total * total))

    # підінтегральний вираз y·x', причому x' = y, тому це квадрат a − h·cos t
    def _integrand_y_dx(self, t):
        return self.point(t)[1] * self.d1(t)[0]

    # площа під однією аркою чисельно, t від 0 до 2π
    def arch_area_numeric(self):
        return simpson(self._integrand_y_dx, 0.0, TWO_PI, 2000)

    # площа під однією аркою точно: A = 2πa² + πh²
    def arch_area_exact(self):
        return math.pi * (2.0 * self.a ** 2 + self.h ** 2)

    # параметр самоперетину T: a·T = h·sin T, бісекція на проміжку від 0 до π
    def loop_half_parameter(self):
        a, h = self.a, self.h
        lo, hi = 1e-12, math.pi
        for _ in range(100):
            mid = (lo + hi) / 2.0
            if a * mid - h * math.sin(mid) < 0:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2.0

    # площа петлі чисельно, інтеграл від −T до T
    def loop_area_numeric(self):
        if not self.is_prolate():
            return None
        half = self.loop_half_parameter()
        return simpson(self._integrand_y_dx, -half, half, 2000)

    # площа петлі точно: A = 2a²T − 4ah·sin T + h²·T + h²·sin 2T / 2
    def loop_area_exact(self):
        if not self.is_prolate():
            return None
        a, h = self.a, self.h
        half = self.loop_half_parameter()
        return (2.0 * a * a * half - 4.0 * a * h * math.sin(half)
                + h * h * (half + math.sin(2.0 * half) / 2.0))


# ---------------- евклідові перетворення ----------------

# рух площини в однорідних координатах: поворот на θ і зсув tx, ty
class RigidMotion:
    def __init__(self):
        self.reset()

    def reset(self):
        self.theta = 0.0
        self.tx = 0.0
        self.ty = 0.0

    # x' = x·cos θ − y·sin θ + tx, y' = x·sin θ + y·cos θ + ty
    def apply(self, x, y):
        c, s = math.cos(self.theta), math.sin(self.theta)
        return c * x - s * y + self.tx, s * x + c * y + self.ty

    # для вектора лише поворот, без зсуву
    def rotate_vector(self, vx, vy):
        c, s = math.cos(self.theta), math.sin(self.theta)
        return c * vx - s * vy, s * vx + c * vy

    # множення на матрицю зсуву додає Δx і Δy
    def translate(self, dx, dy):
        self.tx += dx
        self.ty += dy

    # поворот навколо cx, cy: перенесення центра в початок, поворот, перенесення назад
    def rotate_about(self, cx, cy, alpha):
        c, s = math.cos(alpha), math.sin(alpha)
        dx, dy = self.tx - cx, self.ty - cy
        self.tx = cx + c * dx - s * dy
        self.ty = cy + s * dx + c * dy
        self.theta += alpha


# ---------------- анімація параметра ----------------

# параметр, що плавно змінюється між low і high через синус фази
class Pulse:
    def __init__(self, target, enabled, low, high, speed, minimum):
        self.target = target
        self.enabled = enabled
        self.low = low
        self.high = high
        self.speed = speed
        self.minimum = minimum
        self.phase = 0.0


# ---------------- камера ----------------

# масштаб і зсув камери, перехід світ ↔ екран
class Viewport:
    def __init__(self):
        self.scale = 1.9
        self.offx = 0.0
        self.offy = 0.0
        self.cw = 800
        self.ch = 600

    # y віднімається, бо у canvas вісь Y спрямована вниз
    def to_screen(self, x, y):
        return (self.cw / 2 + (x - self.offx) * self.scale,
                self.ch / 2 - (y - self.offy) * self.scale)

    def to_world(self, sx, sy):
        return ((sx - self.cw / 2) / self.scale + self.offx,
                -(sy - self.ch / 2) / self.scale + self.offy)


# ---------------- програма ----------------

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("ЛР2 — Трохоїда, варіант 16")
        try:
            self.state("zoomed")
        except tk.TclError:
            self.geometry("1400x850")

        self.curve = Trochoid()
        self.motion = RigidMotion()
        self.view = Viewport()
        self.t_min, self.t_max = DEFAULT_T_MIN, DEFAULT_T_MAX

        # побудова
        self.show_grid = tk.BooleanVar(value=True)
        self.show_axes = tk.BooleanVar(value=True)
        self.grid_step = tk.DoubleVar(value=50.0)
        self.curve_color = tk.StringVar(value="black")
        self.curve_width = tk.DoubleVar(value=2.0)

        # параметри кривої
        self.a_var = tk.DoubleVar(value=DEFAULT_A)
        self.h_var = tk.DoubleVar(value=DEFAULT_H)
        self.tmin_var = tk.DoubleVar(value=round(DEFAULT_T_MIN, 4))
        self.tmax_var = tk.DoubleVar(value=round(DEFAULT_T_MAX, 4))
        self.dt_var = tk.DoubleVar(value=DEFAULT_DT)
        self.type_text = tk.StringVar()

        # анімація параметрів h і a
        self.h_anim = tk.BooleanVar(value=False)
        self.h_low = tk.DoubleVar(value=10.0)
        self.h_high = tk.DoubleVar(value=70.0)
        self.h_speed = tk.DoubleVar(value=0.02)
        self.a_anim = tk.BooleanVar(value=False)
        self.a_low = tk.DoubleVar(value=30.0)
        self.a_high = tk.DoubleVar(value=50.0)
        self.a_speed = tk.DoubleVar(value=0.015)
        self.pulses = [Pulse(self.h_var, self.h_anim, self.h_low, self.h_high, self.h_speed, 0.0),
                       Pulse(self.a_var, self.a_anim, self.a_low, self.a_high, self.a_speed, 1.0)]

        # кочення кола
        self.roll_dir = 0
        self.roll_speed = tk.DoubleVar(value=0.02)
        self.show_circle = tk.BooleanVar(value=True)

        # дотична та нормаль
        self.t0_var = tk.DoubleVar(value=DEFAULT_T0)
        self.show_tangent = tk.BooleanVar(value=True)
        self.show_normal = tk.BooleanVar(value=True)
        self.seg_len = tk.DoubleVar(value=160.0)
        self.point_text = tk.StringVar()

        # розрахунки
        self.show_inflections = tk.BooleanVar(value=True)

        # евклідові перетворення
        self.shift_dx = tk.DoubleVar(value=50.0)
        self.shift_dy = tk.DoubleVar(value=0.0)
        self.rot_cx = tk.DoubleVar(value=100.0)
        self.rot_cy = tk.DoubleVar(value=-50.0)
        self.show_center = tk.BooleanVar(value=True)
        self.rot_angle = tk.DoubleVar(value=30.0)
        self.state_text = tk.StringVar()

        # стан миші та перемальовування
        self._drag = None
        self._last_mouse = (0, 0)
        self._redraw_pending = False
        self._range_set = False
        self._wrap_labels = []

        self._build_ui()
        self._bind_canvas()
        for var in (self.show_grid, self.show_axes, self.grid_step, self.curve_color, self.curve_width,
                    self.a_var, self.h_var, self.tmin_var, self.tmax_var, self.dt_var, self.show_circle,
                    self.t0_var, self.show_tangent, self.show_normal, self.seg_len, self.show_inflections,
                    self.rot_cx, self.rot_cy, self.show_center):
            var.trace_add("write", self._schedule_redraw)
        self._schedule_redraw()
        self.after(16, self._tick)

    # ---------------- допоміжні методи ----------------

    # безпечне читання числа з поля, поки користувач друкує
    @staticmethod
    def _num(var, default):
        try:
            v = float(var.get())
        except (tk.TclError, ValueError):
            return default
        return v if math.isfinite(v) else default

    # підпис з автоматичним перенесенням рядків
    def _label(self, parent, textvariable=None, text=""):
        lbl = ttk.Label(parent, text=text, textvariable=textvariable, justify="left")
        lbl.pack(anchor="w", fill="x", padx=6, pady=3)
        self._wrap_labels.append(lbl)
        return lbl

    @staticmethod
    def _row(parent):
        r = ttk.Frame(parent)
        r.pack(fill=tk.X, padx=6, pady=3)
        return r

    @staticmethod
    def _spin(parent, var, lo, hi, inc, width=7):
        sp = ttk.Spinbox(parent, from_=lo, to=hi, increment=inc, textvariable=var, width=width)
        sp.pack(side="left", padx=4)
        return sp

    # округлення значення повзунка, запис лише якщо воно справді змінилось
    def _round_into(self, var, value):
        value = round(float(value), 2)
        if self._num(var, None) != value:
            var.set(value)

    # повзунок разом із полем для точного значення
    def _slider(self, parent, title, var, lo, hi, inc):
        r = self._row(parent)
        ttk.Label(r, text=title, width=3).pack(side="left")
        scale = ttk.Scale(r, from_=lo, to=hi, variable=var, orient="horizontal",
                          command=lambda v: self._round_into(var, v))
        scale.pack(side="left", fill="x", expand=True, padx=4)
        self._spin(r, var, lo, hi, inc, 8)
        return scale

    # ---------------- побудова інтерфейсу ----------------

    def _build_ui(self):
        paned = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True)

        self.canvas = tk.Canvas(paned, bg="white", highlightthickness=0)
        paned.add(self.canvas, weight=4)

        # права панель з прокруткою
        right = ttk.Frame(paned, width=500)
        paned.add(right, weight=1)
        sc = tk.Canvas(right, highlightthickness=0)
        vs = ttk.Scrollbar(right, orient="vertical", command=sc.yview)
        wrap = ttk.Frame(sc)
        wrap.bind("<Configure>", lambda e: sc.configure(scrollregion=sc.bbox("all")))
        inner = sc.create_window((0, 0), window=wrap, anchor="nw")

        def sync_width(e):
            sc.itemconfigure(inner, width=e.width)
            for lbl in self._wrap_labels:
                lbl.configure(wraplength=max(50, e.width - 40))

        sc.bind("<Configure>", sync_width)
        sc.configure(yscrollcommand=vs.set)
        sc.pack(side="left", fill="both", expand=True)
        vs.pack(side="right", fill="y")

        # прокрутка панелі колесом, лише коли курсор над нею
        def wheel(e):
            sc.yview_scroll(int(-e.delta / 120), "units")

        right.bind("<Enter>", lambda e: (sc.bind_all("<MouseWheel>", wheel),
                                         sc.bind_all("<Button-4>", lambda ev: sc.yview_scroll(-1, "units")),
                                         sc.bind_all("<Button-5>", lambda ev: sc.yview_scroll(1, "units"))))
        right.bind("<Leave>", lambda e: (sc.unbind_all("<MouseWheel>"), sc.unbind_all("<Button-4>"),
                                         sc.unbind_all("<Button-5>")))

        # --- побудова: осі, сітка, крива ---
        f = ttk.LabelFrame(wrap, text="Побудова: осі, сітка, крива")
        f.pack(fill=tk.X, padx=8, pady=6)
        r = self._row(f)
        ttk.Checkbutton(r, text="Сітка", variable=self.show_grid).pack(side="left")
        ttk.Checkbutton(r, text="Осі", variable=self.show_axes).pack(side="left", padx=8)
        ttk.Label(r, text="Крок сітки:").pack(side="left")
        self._spin(r, self.grid_step, 5, 500, 5)
        r = self._row(f)
        ttk.Label(r, text="Колір:").pack(side="left")
        ttk.Combobox(r, values=COLORS, textvariable=self.curve_color, width=10,
                     state="readonly").pack(side="left", padx=4)
        ttk.Label(r, text="Товщина:").pack(side="left", padx=(8, 0))
        self._spin(r, self.curve_width, 1, 6, 0.5, 5)

        # --- параметри кривої ---
        f = ttk.LabelFrame(wrap, text="Параметри кривої")
        f.pack(fill=tk.X, padx=8, pady=6)
        ttk.Label(f, text="x(t) = a·t − h·sin t,   y(t) = a − h·cos t").pack(anchor="w", padx=6, pady=3)
        self._slider(f, "a:", self.a_var, 1, 200, 1)
        self._slider(f, "h:", self.h_var, 0, 300, 1)
        r = self._row(f)
        ttk.Label(r, text="t від:").pack(side="left")
        self._spin(r, self.tmin_var, -100, 0, 0.5, 9)
        ttk.Label(r, text="до:").pack(side="left", padx=(8, 0))
        self._spin(r, self.tmax_var, 0, 100, 0.5, 9)
        r = self._row(f)
        ttk.Label(r, text="Крок дискретизації Δt:").pack(side="left")
        self._spin(r, self.dt_var, 0.001, 0.2, 0.005, 7)
        self._label(f, textvariable=self.type_text)
        ttk.Button(f, text="Параметри за замовчуванням", command=self._defaults).pack(anchor="w", padx=6, pady=4)

        # --- анімація ---
        f = ttk.LabelFrame(wrap, text="Анімація")
        f.pack(fill=tk.X, padx=8, pady=6)
        for name, pulse in (("h", self.pulses[0]), ("a", self.pulses[1])):
            ttk.Checkbutton(f, text=f"Автозміна параметра {name}", variable=pulse.enabled,
                            command=lambda p=pulse: self._align_pulse(p)).pack(anchor="w", padx=6, pady=(4, 0))
            r = self._row(f)
            ttk.Label(r, text="min:").pack(side="left")
            self._spin(r, pulse.low, 0, 300, 1, 6)
            ttk.Label(r, text="max:").pack(side="left")
            self._spin(r, pulse.high, 0, 300, 1, 6)
            ttk.Label(r, text="швидк.:").pack(side="left")
            self._spin(r, pulse.speed, 0.001, 0.2, 0.005, 6)
        ttk.Separator(f).pack(fill=tk.X, padx=6, pady=6)
        r = self._row(f)
        ttk.Label(r, text="Кочення кола:").pack(side="left")
        ttk.Button(r, text="◀", width=3, command=lambda: self._roll(-1)).pack(side="left", padx=4)
        ttk.Button(r, text="Пауза", command=lambda: self._roll(0)).pack(side="left")
        ttk.Button(r, text="▶", width=3, command=lambda: self._roll(1)).pack(side="left", padx=4)
        self._spin(r, self.roll_speed, 0.005, 0.2, 0.005, 6)
        ttk.Checkbutton(f, text="Показувати коло, що котиться", variable=self.show_circle).pack(anchor="w", padx=6,
                                                                                               pady=3)

        # --- дотична та нормаль ---
        f = ttk.LabelFrame(wrap, text="Дотична та нормаль у точці t0")
        f.pack(fill=tk.X, padx=8, pady=6)
        r = self._row(f)
        ttk.Label(r, text="t0:", width=3).pack(side="left")
        self.t0_scale = ttk.Scale(r, from_=self.t_min, to=self.t_max, variable=self.t0_var, orient="horizontal",
                                  command=lambda v: self._round_into(self.t0_var, v))
        self.t0_scale.pack(side="left", fill="x", expand=True, padx=4)
        self._spin(r, self.t0_var, -100, 100, 0.1, 8)
        r = self._row(f)
        ttk.Checkbutton(r, text="Дотична", variable=self.show_tangent).pack(side="left")
        ttk.Checkbutton(r, text="Нормаль", variable=self.show_normal).pack(side="left", padx=8)
        ttk.Label(r, text="Довжина:").pack(side="left")
        self._spin(r, self.seg_len, 20, 800, 10)
        self._label(f, textvariable=self.point_text)

        # --- асимптоти ---
        f = ttk.LabelFrame(wrap, text="Асимптоти")
        f.pack(fill=tk.X, padx=8, pady=6)
        self._label(f, text="Асимптот немає. Коефіцієнт k = lim y/x = 0, а lim y при t → ±∞ не існує: "
                            "y = a − h·cos t коливається в межах [a − h; a + h]. Лише при h = 0 крива "
                            "вироджується в пряму y = a.")

        # --- розрахунки ---
        f = ttk.LabelFrame(wrap, text="Розрахунки")
        f.pack(fill=tk.X, padx=8, pady=6)
        self.lbl_len = self._label(f)
        self.lbl_area = self._label(f)
        self.lbl_loop = self._label(f)
        self.lbl_infl = self._label(f)
        ttk.Checkbutton(f, text="Позначати точки перегину на кривій",
                        variable=self.show_inflections).pack(anchor="w", padx=6, pady=3)
        self.lbl_rho = self._label(f)

        # --- евклідові перетворення ---
        f = ttk.LabelFrame(wrap, text="Евклідові перетворення")
        f.pack(fill=tk.X, padx=8, pady=6)
        r = self._row(f)
        ttk.Label(r, text="Δx:").pack(side="left")
        self._spin(r, self.shift_dx, -1000, 1000, 10)
        ttk.Label(r, text="Δy:").pack(side="left")
        self._spin(r, self.shift_dy, -1000, 1000, 10)
        ttk.Button(r, text="Зсунути", command=self._translate).pack(side="left", padx=6)
        r = self._row(f)
        ttk.Label(r, text="cx:").pack(side="left")
        self._spin(r, self.rot_cx, -1000, 1000, 10)
        ttk.Label(r, text="cy:").pack(side="left")
        self._spin(r, self.rot_cy, -1000, 1000, 10)
        ttk.Checkbutton(r, text="Показати", variable=self.show_center).pack(side="left", padx=6)
        r = self._row(f)
        ttk.Label(r, text="Кут Δθ, °:").pack(side="left")
        self._spin(r, self.rot_angle, -360, 360, 5)
        ttk.Button(r, text="Повернути", command=self._rotate).pack(side="left", padx=6)
        self._label(f, textvariable=self.state_text)
        ttk.Button(f, text="Скинути перетворення", command=self._reset_motion).pack(anchor="w", padx=6, pady=4)


    # ---------------- події полотна ----------------

    def _bind_canvas(self):
        c = self.canvas
        c.bind("<Configure>", self._schedule_redraw)
        c.bind("<MouseWheel>", lambda e: self._zoom_at(e.x, e.y, 1.1 if e.delta > 0 else 1 / 1.1))
        c.bind("<Button-4>", lambda e: self._zoom_at(e.x, e.y, 1.1))
        c.bind("<Button-5>", lambda e: self._zoom_at(e.x, e.y, 1 / 1.1))
        c.bind("<ButtonPress-1>", self._on_press)
        c.bind("<B1-Motion>", self._on_drag)
        c.bind("<ButtonRelease-1>", self._on_release)
        c.bind("<Button-3>", self._pick_t0)
        c.bind("<Button-2>", self._pick_t0)

    # масштаб із фіксацією точки під курсором
    def _zoom_at(self, sx, sy, k):
        wx, wy = self.view.to_world(sx, sy)
        self.view.scale = max(0.05, min(40.0, self.view.scale * k))
        wx2, wy2 = self.view.to_world(sx, sy)
        self.view.offx += wx - wx2
        self.view.offy += wy - wy2
        self._schedule_redraw()

    # бірюзова точка — образ початку координат кривої
    def _move_handle(self):
        return self.motion.apply(0.0, 0.0)

    def _on_press(self, e):
        mx, my = self.view.to_screen(*self._move_handle())
        rx, ry = self.view.to_screen(self._num(self.rot_cx, 0.0), self._num(self.rot_cy, 0.0))
        if (e.x - mx) ** 2 + (e.y - my) ** 2 <= 100:
            self._drag = "move"
        elif self.show_center.get() and (e.x - rx) ** 2 + (e.y - ry) ** 2 <= 100:
            self._drag = "center"
        else:
            self._drag = "pan"
        self._last_mouse = (e.x, e.y)

    def _on_drag(self, e):
        if self._drag == "pan":
            self.view.offx -= (e.x - self._last_mouse[0]) / self.view.scale
            self.view.offy += (e.y - self._last_mouse[1]) / self.view.scale
            self._last_mouse = (e.x, e.y)
            self._schedule_redraw()
            return
        wx, wy = self.view.to_world(e.x, e.y)
        if self._drag == "move":
            hx, hy = self._move_handle()
            self.motion.translate(wx - hx, wy - hy)
            self._schedule_redraw()
        elif self._drag == "center":
            # змінюється лише центр, крива стоїть на місці
            self.rot_cx.set(round(wx, 2))
            self.rot_cy.set(round(wy, 2))

    def _on_release(self, _):
        self._drag = None

    # вибір t0 як найближчої до кліку точки кривої
    def _pick_t0(self, e):
        best, best_t = math.inf, None
        for t in self._parameter_values(self.t_min, self.t_max):
            sx, sy = self.view.to_screen(*self.motion.apply(*self.curve.point(t)))
            d = (sx - e.x) ** 2 + (sy - e.y) ** 2
            if d < best:
                best, best_t = d, t
        if best_t is not None:
            self.roll_dir = 0
            self.t0_var.set(round(best_t, 4))

    # ---------------- кнопки ----------------

    # повернення стандартних параметрів кривої та анімації
    def _defaults(self):
        self.h_anim.set(False)
        self.a_anim.set(False)
        self.roll_dir = 0
        self.a_var.set(DEFAULT_A)
        self.h_var.set(DEFAULT_H)
        self.tmin_var.set(round(DEFAULT_T_MIN, 4))
        self.tmax_var.set(round(DEFAULT_T_MAX, 4))
        self.dt_var.set(DEFAULT_DT)
        self.t0_var.set(DEFAULT_T0)

    # напрям кочення: −1 назад, 0 пауза, 1 вперед
    def _roll(self, direction):
        self.roll_dir = direction

    def _translate(self):
        self.motion.translate(self._num(self.shift_dx, 0.0), self._num(self.shift_dy, 0.0))
        self._schedule_redraw()

    def _rotate(self):
        alpha = math.radians(self._num(self.rot_angle, 0.0))
        self.motion.rotate_about(self._num(self.rot_cx, 0.0), self._num(self.rot_cy, 0.0), alpha)
        self._schedule_redraw()

    def _reset_motion(self):
        self.motion.reset()
        self._schedule_redraw()

    # ---------------- таймер анімації ----------------

    # межі анімації, впорядковані та не менші за допустимий мінімум параметра
    def _pulse_bounds(self, pulse):
        low = max(pulse.minimum, self._num(pulse.low, pulse.minimum))
        high = max(pulse.minimum, self._num(pulse.high, pulse.minimum))
        return min(low, high), max(low, high)

    # початкова фаза через арксинус, щоб анімація стартувала від поточного значення
    def _align_pulse(self, pulse):
        low, high = self._pulse_bounds(pulse)
        amp = (high - low) / 2.0
        if amp <= 0:
            pulse.phase = 0.0
            return
        r = (self._num(pulse.target, low) - (low + high) / 2.0) / amp
        pulse.phase = math.asin(max(-1.0, min(1.0, r)))

    # кожні 16 мс фаза зростає, параметр рахується через синус, а t0 рухається при коченні
    def _tick(self):
        for pulse in self.pulses:
            if not pulse.enabled.get():
                continue
            low, high = self._pulse_bounds(pulse)
            pulse.phase += self._num(pulse.speed, 0.0)
            value = (low + high) / 2.0 + (high - low) / 2.0 * math.sin(pulse.phase)
            pulse.target.set(round(value, 2))
        if self.roll_dir:
            t0 = self._num(self.t0_var, DEFAULT_T0) + self.roll_dir * self._num(self.roll_speed, 0.02)
            if t0 > self.t_max:
                t0 = self.t_min
            elif t0 < self.t_min:
                t0 = self.t_max
            self.t0_var.set(round(t0, 4))
        self.after(16, self._tick)

    # ---------------- від формул до екрана ----------------

    # читання полів у модель із перевіркою значень
    def _read_params(self):
        self.curve.a = max(1.0, self._num(self.a_var, DEFAULT_A))
        self.curve.h = max(0.0, self._num(self.h_var, DEFAULT_H))
        t1 = self._num(self.tmin_var, DEFAULT_T_MIN)
        t2 = self._num(self.tmax_var, DEFAULT_T_MAX)
        if t2 <= t1:
            t1, t2 = DEFAULT_T_MIN, DEFAULT_T_MAX
        # межі повзунка змінюються лише тоді, коли змінився діапазон, інакше перемальовування зациклюється
        if (t1, t2) != (self.t_min, self.t_max) or not self._range_set:
            self._range_set = True
            self.t0_scale.configure(from_=t1, to=t2)
        self.t_min, self.t_max = t1, t2
        self.view.cw = self.canvas.winfo_width()
        self.view.ch = self.canvas.winfo_height()

    # поточне t0 в межах діапазону
    def _current_t0(self):
        return max(self.t_min, min(self.t_max, self._num(self.t0_var, DEFAULT_T0)))

    # рівні кроки Δt, але не більше 8000 точок
    def _parameter_values(self, t1, t2):
        dt = max(1e-4, self._num(self.dt_var, DEFAULT_DT))
        n = int(min(MAX_POINTS, max(2, math.ceil((t2 - t1) / dt))))
        step = (t2 - t1) / n
        return [t1 + i * step for i in range(n + 1)]

    # формула кривої, потім евклідове перетворення, потім пікселі
    def _curve_to_screen(self, ts):
        coords = []
        for t in ts:
            x, y = self.motion.apply(*self.curve.point(t))
            coords.extend(self.view.to_screen(x, y))
        return coords

    # кілька змін поспіль дають одне перемальовування
    def _schedule_redraw(self, *_):
        if not self._redraw_pending:
            self._redraw_pending = True
            self.after_idle(self._redraw)

    # шари по черзі: сітка, осі, крива, коло, точки перегину, дотична й нормаль, ручки
    def _redraw(self):
        self._redraw_pending = False
        self._read_params()
        self.canvas.delete("all")
        if self.show_grid.get():
            self._draw_grid()
        if self.show_axes.get():
            self._draw_axes()
        self._draw_curve()
        self._draw_rolling_circle()
        self._draw_inflections()
        self._draw_tangent_normal()
        self._draw_handles()
        self._update_info()

    # ---------------- шари малювання ----------------

    def _draw_grid(self):
        c = self.canvas
        step = self._num(self.grid_step, 50.0)
        if step <= 0 or step * self.view.scale < 4:
            return
        x0, y0 = self.view.to_world(0, self.view.ch)
        x1, y1 = self.view.to_world(self.view.cw, 0)
        k = math.floor(x0 / step)
        while k * step <= x1:
            sx, _ = self.view.to_screen(k * step, 0)
            c.create_line(sx, 0, sx, self.view.ch, fill="#e9ecef")
            k += 1
        k = math.floor(y0 / step)
        while k * step <= y1:
            _, sy = self.view.to_screen(0, k * step)
            c.create_line(0, sy, self.view.cw, sy, fill="#e9ecef")
            k += 1

    # осі зі стрілками та підписами поділок
    def _draw_axes(self):
        c = self.canvas
        cw, ch = self.view.cw, self.view.ch
        ox, oy = self.view.to_screen(0, 0)
        oy = max(10, min(ch - 10, oy))
        ox = max(10, min(cw - 10, ox))
        c.create_line(0, oy, cw - 2, oy, fill="#d62828", width=2, arrow="last", arrowshape=(12, 14, 6))
        c.create_text(cw - 14, oy - 16, text="X", fill="#d62828", font=("Segoe UI", 11, "bold"))
        c.create_line(ox, ch, ox, 2, fill="#2a9d8f", width=2, arrow="last", arrowshape=(12, 14, 6))
        c.create_text(ox + 16, 16, text="Y", fill="#2a9d8f", font=("Segoe UI", 11, "bold"))
        step = self._num(self.grid_step, 50.0)
        if step <= 0 or step * self.view.scale < 30:
            return
        x0, y0 = self.view.to_world(0, ch)
        x1, y1 = self.view.to_world(cw, 0)
        k = math.floor(x0 / step)
        while k * step <= x1:
            if k:
                sx, _ = self.view.to_screen(k * step, 0)
                c.create_text(sx, oy + 14, text=f"{k * step:g}", fill="#868e96", font=("Segoe UI", 8))
            k += 1
        k = math.floor(y0 / step)
        while k * step <= y1:
            if k:
                _, sy = self.view.to_screen(0, k * step)
                c.create_text(ox - 8, sy, text=f"{k * step:g}", fill="#868e96", anchor="e", font=("Segoe UI", 8))
            k += 1

    def _draw_curve(self):
        coords = self._curve_to_screen(self._parameter_values(self.t_min, self.t_max))
        if len(coords) >= 4:
            self.canvas.create_line(*coords, fill=self.curve_color.get(), width=self._num(self.curve_width, 2.0))

    # коло з центром a·t0, a, спиця до точки і синій слід уже описаної частини
    def _draw_rolling_circle(self):
        if not self.show_circle.get():
            return
        c = self.canvas
        a, t0 = self.curve.a, self._current_t0()
        # пряма, по якій котиться колесо
        p1 = self.view.to_screen(*self.motion.apply(a * self.t_min - a, 0.0))
        p2 = self.view.to_screen(*self.motion.apply(a * self.t_max + a, 0.0))
        c.create_line(*p1, *p2, fill="#adb5bd", dash=(5, 4))
        if t0 > self.t_min:
            coords = self._curve_to_screen(self._parameter_values(self.t_min, t0))
            c.create_line(*coords, fill="#1c7ed6", width=self._num(self.curve_width, 2.0))
        cx, cy = self.view.to_screen(*self.motion.apply(a * t0, a))
        r = a * self.view.scale
        c.create_oval(cx - r, cy - r, cx + r, cy + r, outline="#6c757d")
        px, py = self.view.to_screen(*self.motion.apply(*self.curve.point(t0)))
        c.create_line(cx, cy, px, py, fill="#6c757d")
        c.create_oval(cx - 2.5, cy - 2.5, cx + 2.5, cy + 2.5, fill="#6c757d", outline="")

    # помаранчеві ромби в точках перегину
    def _draw_inflections(self):
        if not self.show_inflections.get():
            return
        for t in self.curve.inflections(self.t_min, self.t_max):
            sx, sy = self.view.to_screen(*self.motion.apply(*self.curve.point(t)))
            self.canvas.create_polygon(sx, sy - 7, sx + 7, sy, sx, sy + 7, sx - 7, sy,
                                       fill="#fd7e14", outline="#a04800")

    # дотична — одиничний вектор швидкості, нормаль — він же, повернутий на 90°
    def _draw_tangent_normal(self):
        t0 = self._current_t0()
        px, py = self.motion.apply(*self.curve.point(t0))
        show_tangent = self.show_tangent.get()
        show_normal = self.show_normal.get()
        half = self._num(self.seg_len, 160.0) / 2.0
        ux, uy = self.motion.rotate_vector(*self.curve.unit_tangent(t0))
        if show_tangent:
            self.canvas.create_line(*self.view.to_screen(px - ux * half, py - uy * half),
                                    *self.view.to_screen(px + ux * half, py + uy * half),
                                    fill="#dc3545", width=2)
        if show_normal:
            self.canvas.create_line(*self.view.to_screen(px + uy * half, py - ux * half),
                                    *self.view.to_screen(px - uy * half, py + ux * half),
                                    fill="#198754", width=2)
        sx, sy = self.view.to_screen(px, py)
        self.canvas.create_oval(sx - 6, sy - 6, sx + 6, sy + 6, fill="#1c7ed6", outline="")
        self.point_text.set(f"t0 = {t0:.4f};  P = ({px:.3f}; {py:.3f})")

    # бірюзова точка зсуву та жовта точка центру обертання
    def _draw_handles(self):
        c = self.canvas
        mx, my = self.view.to_screen(*self._move_handle())
        c.create_oval(mx - 7, my - 7, mx + 7, my + 7, fill="#22c3d6", outline="#0e7f90")
        if self.show_center.get():
            rx, ry = self.view.to_screen(self._num(self.rot_cx, 0.0), self._num(self.rot_cy, 0.0))
            c.create_oval(rx - 7, ry - 7, rx + 7, ry + 7, fill="#fcc419", outline="#b38f00")

    # ---------------- вивід розрахунків ----------------

    def _update_info(self):
        cv = self.curve
        self.type_text.set(f"Тип: {cv.kind()}")

        length = cv.arc_length(self.t_min, self.t_max)
        self.lbl_len.configure(
            text=f"Довжина дуги на [{self.t_min:.3f}; {self.t_max:.3f}]: L ≈ {length:.3f}\n"
                 f"За період: чисельно {cv.arc_length(0.0, TWO_PI):.3f}; "
                 f"4(a+h)·E(m) = {cv.period_length_exact():.3f}")

        area_text = (f"Площа під аркою, t ∈ [0; 2π]: чисельно {cv.arch_area_numeric():.3f}; "
                     f"π(2a²+h²) = {cv.arch_area_exact():.3f}")
        if cv.is_prolate():
            area_text += "\nПри h > a інтеграл ∫y dx враховує знак і кратність обходу."
        self.lbl_area.configure(text=area_text)

        if cv.is_prolate():
            self.lbl_loop.configure(text=f"Площа однієї петлі: чисельно {cv.loop_area_numeric():.3f}; "
                                         f"формула {cv.loop_area_exact():.3f}")
        else:
            self.lbl_loop.configure(text="Площа петлі: петель немає, h ≤ a.")

        if cv.is_curtate():
            points = cv.inflections(self.t_min, self.t_max)
            self.lbl_infl.configure(
                text=f"Точки перегину: cos t = h/a, t = ±{math.acos(cv.h / cv.a):.4f} + 2πk, "
                     f"y = a − h²/a = {cv.a - cv.h * cv.h / cv.a:.3f}. У діапазоні t: {len(points)} шт.")
        elif cv.is_cycloid():
            self.lbl_infl.configure(text="Точки перегину: немає, при h = a у точках t = 2πk вістря, а не перегин.")
        elif cv.h == 0:
            self.lbl_infl.configure(text="Точки перегину: немає, при h = 0 крива — пряма.")
        else:
            self.lbl_infl.configure(text="Точки перегину: немає, бо при h > a рівняння a·cos t = h "
                                         "не має розв'язків.")

        t0 = self._current_t0()
        rho = cv.radius_of_curvature(t0)
        if cv.is_cusp(t0):
            rho_text = "0 (вістря)"
        elif math.isinf(rho):
            rho_text = "∞ (кривина нульова)"
        else:
            rho_text = f"{rho:.3f}"
        self.lbl_rho.configure(text=f"Радіус кривини: ρ(t0) = {rho_text}")

        self.state_text.set(f"Поточний стан: θ = {math.degrees(self.motion.theta):.2f}°, "
                            f"T = ({self.motion.tx:.2f}; {self.motion.ty:.2f})")


if __name__ == "__main__":
    App().mainloop()
