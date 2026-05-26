f_data = [10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, 60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, 85.0, 87.5, 90.0, 92.5, 95.0, 97.5, 100.0, 102.5, 105.0, 107.5]
V_data = [0.842, 0.911, 0.986, 1.062, 1.143, 1.227, 1.314, 1.401, 1.482, 1.551, 1.216, 1.048, 0.866, 0.689, 0.521, 0.364, 0.223, 0.103, 0.012, -0.041, -0.057, -0.034, 0.018, 0.096, 0.197, 0.318, 0.452, 0.579, 0.700, 0.809, 0.611, 0.688, 0.756, 0.811, 0.856, 0.894, 0.926, 0.954, 0.980, 1.004]
h = 2.5

# =================================================================
# 1. INTERVALOS DONDE V(f) CAMBIA DE SIGNO
# =================================================================
print("1. INTERVALOS DONDE V(f) CAMBIA DE SIGNO")
print("-----------------------------------------------------------------")
intv = []
for i in range(len(V_data) - 1):
    if V_data[i] * V_data[i+1] < 0:
        intv.append((f_data[i], f_data[i+1], V_data[i], V_data[i+1]))

print(f"Se detectaron {len(intv)} cruces por cero en la tabla:")
print(f"  > Primer Intervalo:  [{intv[0][0]} kHz, {intv[0][1]} kHz]  (V pasa de positivo a negativo)")
print(f"  > Segundo Intervalo: [{intv[1][0]} kHz, {intv[1][1]} kHz]  (V pasa de negativo a positivo)")
print("\n" + "="*65 + "\n")


# =================================================================
# 2. BISECCIÓN PARA LOCALIZAR LAS RAÍCES (DATOS DE LA TABLA)
# =================================================================
print("2. BISECCIÓN PARA LOCALIZAR LAS RAÍCES (DATOS DE LA TABLA)")
print("-----------------------------------------------------------------")

# Leyenda explicativa global al inicio del método
print("Nomenclatura del método de bisección para todas las tablas:")
print("  a    = Límite inferior del intervalo de búsqueda.")
print("  b    = Límite superior del intervalo de búsqueda.")
print("  m    = (a + b) / 2, punto medio del intervalo actual.")
print("  V(m) = Voltaje evaluado en el punto medio.")
print("-----------------------------------------------------------------")

def bisec_tabla(a, b, va, vb, num_raiz):
    print(f"Proceso de Bisección - Alarma {num_raiz} en Tabla:")
    print("+-----------+-----------+-----------+-----------+-----------+")
    print("| Intervalo |     a     |     b     |     m     |   V(m)    |")
    print("+-----------+-----------+-----------+-----------+-----------+")
    for i in range(1, 4):
        m = (a + b) / 2.0
        vm = va + ((vb - va) / (b - a)) * (m - a)
        print(f"|     {i}     |  {a:>6.2f}   |  {b:>6.2f}   |  {m:>7.4f}  |  {vm:>7.4f}  |")
        if va * vm < 0: b = m
        else: a, va = m, vm
    print("+-----------+-----------+-----------+-----------+-----------+")
    return m

# Ejecución de tablas consecutivas
r1_tab = bisec_tabla(intv[0][0], intv[0][1], intv[0][2], intv[0][3], "1")
print()
r2_tab = bisec_tabla(intv[1][0], intv[1][1], intv[1][2], intv[1][3], "2")
print("\n" + "="*65 + "\n")


# =================================================================
# 3. REFINANDO LA RAÍZ CON SPLINE CÚBICO
# =================================================================
print("3. REFINANDO LA RAÍZ CON SPLINE CÚBICO NATURAL")
print("-----------------------------------------------------------------")

n = len(f_data) - 1
A, B, C, D = [0.0]*(n+1), [1.0]*(n+1), [0.0]*(n+1), [0.0]*(n+1)
for i in range(1, n):
    A[i], B[i], C[i] = h/6.0, (2*h)/3.0, h/6.0
    D[i] = ((V_data[i+1] - V_data[i])/h) - ((V_data[i] - V_data[i-1])/h)

cp, dp = [0.0]*(n+1), [0.0]*(n+1)
cp[0], dp[0] = C[0]/B[0], D[0]/B[0]
for i in range(1, n+1):
    den = B[i] - A[i]*cp[i-1]
    if i < n: cp[i] = C[i]/den
    dp[i] = (D[i] - A[i]*dp[i-1])/den
M = [0.0]*(n+1)
M[-1] = dp[-1]
for i in range(n-1, -1, -1): M[i] = dp[i] - cp[i]*M[i+1]

def eval_S(x_obj):
    for i in range(n):
        if f_data[i] <= x_obj <= f_data[i+1]:
            t1 = M[i]*((f_data[i+1]-x_obj)**3)/(6*h) + M[i+1]*((x_obj-f_data[i])**3)/(6*h)
            t3 = (V_data[i]/h - M[i]*h/6)*(f_data[i+1]-x_obj) + (V_data[i+1]/h - M[i+1]*h/6)*(x_obj-f_data[i])
            return t1 + t3

def bisec_spline(a, b, num_raiz):
    print(f"Refinamiento con Spline - Alarma {num_raiz}:")
    print("+-----------+-----------+-----------+-----------+-----------+")
    print("| Intervalo |     a     |     b     |     m     |   V(m)    |")
    print("+-----------+-----------+-----------+-----------+-----------+")
    for i in range(1, 5):
        m = (a + b) / 2.0
        vm = eval_S(m)
        print(f"|     {i}     |  {a:>6.2f}   |  {b:>6.2f}   |  {m:>7.4f}  |  {vm:>7.4f}  |")
        if eval_S(a) * vm < 0: b = m
        else: a = m
    print("+-----------+-----------+-----------+-----------+-----------+")
    return m

r1_spl = bisec_spline(intv[0][0], intv[0][1], "1")
print()
r2_spl = bisec_spline(intv[1][0], intv[1][1], "2")

# Tabla Comparativa Final
print("\nTABLA COMPARATIVA DE ALARMAS (RESULTADO REFINADO)")
print("+-------------------+---------------------+---------------------+")
print(f"| Frecuencia Alarma | Bisección en Tabla  |   Refinada Spline   |")
print("+-------------------+---------------------+---------------------+")
print(f"|  1° Cruce por 0   |     {r1_tab:.4f} kHz     |     {r1_spl:.4f} kHz     |")
print(f"|  2° Cruce por 0   |     {r2_tab:.4f} kHz     |     {r2_spl:.4f} kHz     |")
print("+-------------------+---------------------+---------------------+")