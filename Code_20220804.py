import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.special import logsumexp
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

def DiffEq_cf(d_list, t):
    time = t

    j = exp[0, 2]  # the reactor turnover rate
    Cinfluent = exp[0, 3]  # the inflow As(III) concentration
    k = exp[0, 4]  # rate constant
    Cs = exp[0, 5]  # sorbent concentrations
    KF = exp[0, 6]  # Freundlich constant
    n = exp[0, 7]  # experimentally determined constants
    Qmax = exp[0, 8]
    KL = exp[0, 9]
    isotherm = exp[0, 10]
    order = exp[0, 11]
    Figure = exp[0, 12]

    Ct, qt = d_list  # concentration of aqueous As(III) and the quantity of arsenic adsorbed at time t
    if Ct < 0.001:
        Ct = 0.001  # set a minimum concentration of sorbate in the reactor to prevent the Freundlich adsorption isotherm from running an error.

    if Ct >= Cinfluent:
        Ct = Cinfluent  # control incase ode15s time intervals are too large

    Ct_mgL = Ct / 1000
    qt_mgg = qt / (Cs * 1000)
    if isotherm == 1:
        # Freundlich
        rate_ads = 1000 * k * Ct_mgL * Cs * (
                    (1 - (qt_mgg / (KF * (Ct_mgL ** (1 / n))))) ** order)  # batch treatment reactor design

    elif isotherm == 2:
        # Langmuir
        rate_ads = 1000 * k * Ct_mgL * Cs * ((1 - (qt_mgg / ((Qmax * KL * Ct_mgL) / (1 + KL * Ct_mgL)))) ** order)

    rate_influx = j * Cinfluent  # calculate the rate of sorbate influx (continuous-flow systems only)
    rate_outflux = j * Ct  # calculate the rate of sorbate outflux (continuous-flow systems only)
    dCdt = [0, 0]
    dCdt[0] = -rate_ads + rate_influx - rate_outflux  # continuous treatment reactor design

    if Figure == True:

        dCdt[1] = rate_ads / Cs * 1000
    elif Figure == False:

        dCdt[1] = rate_ads

    return dCdt


def DiffEq_batch(conditions, t):
    time = t

    j = exp[0][2]  # the reactor turnover rate
    Cinfluent = exp[0][3]  # the inflow As(III) concentration
    k = exp[0][4]  # rate constant
    Cs1 = exp[0][5] # sorbent concentrations
    KF = exp[0][6]  # Freundlich constant
    n = exp[0][7]  # experimentally determined constants
    Qmax = exp[0][8]
    KL = exp[0][9]
    isotherm = exp[0][10]

    order = exp[0][11]

    Ct = conditions[0]  # concentration of aqueous As(III) at time t
    qt[i] = conditions[1]  # the quantity of arsenic adsorbed at time t
    Ct_mgL = Ct / 1000
    qt_mgg = qt[i] / (Cs1 * 1000)

    if isotherm == 1:
        # Freundlich
        rate_ads = 1000 * k * Ct_mgL * Cs1 * (
                    (1 - (qt_mgg / (KF * (Ct_mgL ** (1 / n))))) ** order)  # batch treatment reactor design

    elif isotherm == 2:
        # Langmuir
        rate_ads = 1000 * k * Ct_mgL * Cs1 * ((1 - (qt_mgg / ((Qmax * KL * Ct_mgL) / (1 + KL * Ct_mgL)))) ** order)

    if Ct < 0.000001:
        Ct = 0.00001
        rate_ads = 0

    rate_influx = j * Cinfluent
    rate_outflux = j * Ct
    dCdt = [-rate_ads + rate_influx - rate_outflux, rate_ads]  # batch treatment reactor design (when j=0)

    qt[i] = qt[i] + rate_ads

    return dCdt


def Result_Processing_batch(C_in, Cs0, days, ord_in=2, isotherms="Freundlich", k_in=0.1111, KF_in=5.10, n_in=2.63,
                      Qmax_in=29.44, KL_in=0.11):
    global i
    global exp
    global qt
    qt = [0, 0, 0, 0, 0, 0, 0]
    number_of_experiments = 7
    number_of_variables = 12
    results_table = [([0] * number_of_experiments) for i in range(10)]
    if isotherms == "Freundlich":
        isoth_in = 1
    elif isotherms == "Langmuir":
        isoth_in = 2

    for day in range(1, days + 1):
        C_init = C_in  # sorbent concentrations (initial value)
        q_init = qt    # the quantity of arsenic adsorbed at time t (initial value)
        Cs = Cs0       # sorbent concentrations
        Cinfluent = 0  # the inflow As(III) concentration
        j = 0          # the reactor turnover rate
        k = k_in       # rate constant

        KF = KF_in     # Freundlich constant
        Qmax = Qmax_in
        isotherm = isoth_in
        order = ord_in
        KL = KL_in
        n = n_in      # experimentally determined constants

        data_collect = [0, 1, 2, 3, 4,
                        5, 6, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 200, 210, 220, 230,
                        240, 270, 300, 330, 360, 390, 420, 450, 480, 500, 570,
                        600, 660, 690, 720, 750, 780, 840, 870, 900, 930, 960, 990, 1000, 1050, 1080, 1110, 1140,
                        1170, 1200, 1260, 1290, 1320, 1350, 1380, 1410, 1440]
        t_end = 1440

        t_steps = len(data_collect)
        t_step = t_end / t_steps
        tt = list(np.arange(0, t_end + 1, t_step))

        results_Ct = np.zeros((number_of_experiments, len(data_collect)))
        results_qt = np.zeros((number_of_experiments, len(data_collect)))

        experiments = [[0] * number_of_variables for _ in range(number_of_experiments)]
        Cs_col = [0] * number_of_experiments
        for i in range(number_of_experiments):
            experiments[i][0] = C_init
            experiments[i][1] = q_init[i]
            experiments[i][2] = j
            experiments[i][3] = Cinfluent
            experiments[i][4] = k
            experiments[i][5] = Cs * 10 ** (i - 4)  # exponentially increasing sorbent concentration
            experiments[i][6] = KF
            experiments[i][7] = n
            experiments[i][8] = Qmax
            experiments[i][9] = KL
            experiments[i][10] = isotherm
            experiments[i][11] = order
            Cs_col[i] = experiments[i][5]

        for i in range(number_of_experiments):

            exp = [[0] * number_of_variables for _ in range(1)]
            exp[0][:] = experiments[i][:]
            exp_C_init = exp[0][0]
            exp_q_init = exp[0][1]

            result = odeint(DiffEq_batch, [exp_C_init, exp_q_init], data_collect, rtol=1e-6)

            results_Ct[i, :] = result[:, 0]

            for j in range(len(results_Ct[i, :])):
                if results_Ct[i, :][j] <= 10 and data_collect[j] == 1:
                    results_table[0][i] += 1

                if results_Ct[i, :][j] <= 10 and data_collect[j] == 2:
                    results_table[1][i] += 1

                if results_Ct[i, :][j] <= 10 and data_collect[j] == 5:
                    results_table[2][i] += 1

                if results_Ct[i, :][j] <= 10 and data_collect[j] in [10]:
                    results_table[3][i] += 1

                if results_Ct[i, :][j] <= 10 and data_collect[j] in [20]:
                    results_table[4][i] += 1

                if results_Ct[i, :][j] <= 10 and data_collect[j] in [50]:
                    results_table[5][i] += 1

                if results_Ct[i, :][j] <= 10 and data_collect[j] in [100]:
                    results_table[6][i] += 1

                if results_Ct[i, :][j] <= 10 and data_collect[j] in [200]:
                    results_table[7][i] += 1

                if results_Ct[i, :][j] <= 10 and data_collect[j] in [500]:
                    results_table[8][i] += 1

                if results_Ct[i, :][j] <= 10 and data_collect[j] in [100]:
                    results_table[9][i] += 1

            results_qt[i, :] = result[:, 1]

    return results_Ct, Cs_col,data_collect, results_table


def Result_Processing_cf(j_flow, C_inf, Cs0, figure=True, ord_in=2, isotherms1="Freundlich",
                      k_in1=0.1111, KF_in1=5.10, n_in1=2.63,
                      Qmax_in1=29.44, KL_in1=0.11):
    number_of_experiments = 1  # each experiment will have a different initial sorbate concentration (C0) 实验次数

    number_of_variables = 13  # counting how many rows in the array we need to store the input parameters for each kinetic plot 参数数量

    experiments = np.zeros((number_of_experiments,
                            number_of_variables))  # each experiment refers to a single kinetic plot  54*8的矩阵（行为第几次试验，列表示参数）

    if isotherms1 == "Freundlich":
        isoth_in = 1
    elif isotherms1 == "Langmuir":
        isoth_in = 2

    # creating the input variables for each kinetic plots

    C_init = 0.01  # ug L-1  - this is the initial concentration of aqueous sorbate in the suspension 这是悬浮液中山梨醇水溶液的初始浓度

    q_init = 0  # ug L-1  - this is the initial concentration of adsorbed sorbate in the suspension 吸收AS的浓度

    Cs = Cs0  # g L-1 - this is the concentration of sorbent  吸附剂浓度

    Cinfluent = C_inf  # ppb or ug L-1 - this is the concentration of sorbate in the influent (for continuous-flow modelling)

    j = j_flow  # this the turn-over frequency, i.e. bed volumes per minute

    k = k_in1  # this is the value of normalised k' (L g-1 min-1)

    KF = KF_in1  # this is the Freundlich constant (mg g-1) used to determine 'qe' at each time step

    n = n_in1 # this is the second parameter for the Freundlich adsorption isotherm (g L-1) used to determine 'qe' at each time step

    Figure = figure

    Qmax = Qmax_in1

    isotherm = isoth_in

    order = ord_in

    KL = KL_in1

    # setting the time intervals upon which data is recorded. This will be
    # overwritten later to avoid wasting computational time after breakthrough
    # has already occurred.

    t_end = 1000000  # 1440 = 1 day
    t_steps = 20000
    t_step = t_end / t_steps
    bv_end = t_end * j
    tt = list(np.arange(0, t_end + 1, t_step))

    results_table = [0] * number_of_experiments

    for i in range(number_of_experiments):
        experiments[i, 0] = C_init
        experiments[i, 1] = q_init
        experiments[i, 2] = j
        experiments[i, 3] = Cinfluent
        experiments[i, 4] = k
        experiments[i, 5] = Cs  # exponentially increasing sorbent concentration
        experiments[i, 6] = KF
        experiments[i, 7] = n
        experiments[i, 8] = Qmax
        experiments[i, 9] = KL
        experiments[i, 10] = isotherm
        experiments[i, 11] = order
        experiments[i, 12] = Figure

    for i in range(number_of_experiments):
        global exp
        exp = np.zeros((1, number_of_variables))

        exp[0, :] = experiments[i, :]
        exp_C_init = exp[0, 0]
        exp_q_init = exp[0, 1]
        t_end = 180000 * exp[0, 5] / (exp[0, 2] * exp[0, 3])  # 50000*Cs/(j*Cinf)
        if t_end < 1000:
            t_end = t_end * 20

        if t_end < 500:
            t_end = t_end * 5

        if j < 0.01:
            t_end = t_end * 10

        if Figure == True:
            stepping = [i for i in range(
                t_steps + 1)]  # collect data at shorter time intervals in the initial stages of the simulation
            stepping[0] = 0
            gradient_1 = 1
            gradient_2 = 30
            for j in range(1, (t_steps + 1)):
                # change the time intervals from being evenly spaced to having a smooth
                # transition from gradient_1 to gradient_2
                stepping[j] = stepping[j - 1] + (gradient_2 * (stepping[j] / t_steps)) + (
                            gradient_1 * ((t_steps - stepping[j]) / t_steps))

            for k in range(1, (t_steps + 1)):
                # normalise the time intervals to 1 and then multiply out by the desired
                # final time
                stepping[k] = (stepping[k] / stepping[-1]) * t_end
            # result = odeint(DiffEq,[exp_C_init,exp_q_init],stepping)
        elif Figure == False:
            stepping = [i for i in range(0, 525600, 1440)]

        result = odeint(DiffEq_cf, [exp_C_init, exp_q_init], stepping, rtol=1e-6)

        counter = 0
        results_Ct = np.zeros((number_of_experiments, len(result[:, 0])))
        results_t = np.zeros((number_of_experiments, len(stepping)))
        results_t[i, :] = stepping
        results_Ct[i, :] = result[:, 0]

        for j in range(len(result[:, 0])):
            if result[:, 0][j] <= 10.5:
                counter += 1

        results_table[i] = [results_t[i, :], results_Ct[i, :]]
    return results_table, counter


def Min_adsorbent_bt(results_counter,Cs_col):
    ads = Cs_col[2:]
    results = [0] * len(ads)
    for h in range(2, len(results_counter)):
        if results_counter[h]==0:
          results_counter[h]=0.0001
        results[h - 2] = (365 / results_counter[h]) * ads[h - 2] * 0.001 * 40

    return min(results)


def Adsorbents_required(C_in1=500, Cs_in=100, j_flow=0.001, ord_in=2, isotherms1="Freundlich",
                        k_in1=0.1111, KF_in1=5.10, n_in1=2.63, Qmax_in1=29.44, KL_in1=0.11):
    number_of_experiments = 7
    Cs_col_cf = [0] * number_of_experiments
    for i in range(number_of_experiments):
        Cs_col_cf[i] = Cs_in * 10 ** (i - 4)

    _, counter_10000 = Result_Processing_cf(j_flow, C_in1, Cs_col_cf[6], False, ord_in, isotherms1,
                                            k_in1, KF_in1, n_in1, Qmax_in1, KL_in1)
    _, counter_1000 = Result_Processing_cf(j_flow, C_in1, Cs_col_cf[5], False, ord_in, isotherms1,
                                           k_in1, KF_in1, n_in1, Qmax_in1, KL_in1)
    _, counter_100 = Result_Processing_cf(j_flow, C_in1, Cs_col_cf[4], False, ord_in, isotherms1,
                                          k_in1, KF_in1, n_in1, Qmax_in1, KL_in1)
    _, counter_10 = Result_Processing_cf(j_flow, C_in1, Cs_col_cf[3], False, ord_in, isotherms1,
                                         k_in1, KF_in1, n_in1, Qmax_in1, KL_in1)
    _, counter_1 = Result_Processing_cf(j_flow, C_in1, Cs_col_cf[2], False, ord_in, isotherms1,
                                        k_in1, KF_in1, n_in1, Qmax_in1, KL_in1)
    ads = Cs_col_cf[2:]

    result_counter = [counter_1, counter_10, counter_100, counter_1000, counter_10000]
    results = [0] * len(result_counter)
    for i in range(len(result_counter)):
        if result_counter[i] == 0:
            result_counter[i] = 0.0001
        results[i] = (365 / result_counter[i]) * ads[i] * 0.001 * 40

    return min(results)


def main_window():
    global var_order
    global Cin
    global Cs
    global k_in

    window = tk.Tk()
    window.title('Predicting Removal of Pollutant from Sorbents')
    window.geometry('450x300+400+200')
    var_order = tk.StringVar()
    var_order.set(1)
    Cin = tk.StringVar()
    Cin.set(500)
    Cs = tk.StringVar()
    Cs.set(100)
    k_in = tk.StringVar()
    k_in.set(0.1111)

    # user information
    l = tk.Label(window, text='Predicting Removal of Arsenic from Sorbents',
                 bg='Blue', font=('Arial', 12), width=40,
                 height=2)

    l.pack()  # Label内容content区域放置位置，自动调节尺寸
    tk.Label(window, text='Initial sorbate/influence concentration ').place(x=1, y= 50)
    tk.Label(window, text='Initial sorbent concentration').place(x=1, y= 80)
    tk.Label(window, text="Rate constant (k')").place(x=1, y=110)
    tk.Label(window, text="Order of pseudo").place(x=1, y=140)


    entry_initial_sorbent_con = tk.Entry(window,textvariable=Cin)
    entry_initial_sorbent_con.place(x=240, y=50)


    entry_initial_sorbate_con = tk.Entry(window,textvariable=Cs)
    entry_initial_sorbate_con.place(x=240, y=80)

    entry_k = tk.Entry(window,textvariable=k_in)
    entry_k.place(x=240, y=110)


    # Define checkboxes
    entry_O = tk.Checkbutton(window, text="pseudo-second order", variable=var_order, onvalue=2, offvalue=1)
    entry_O.place(x=240, y=140)




    # login and sign up button

    btn_bt = tk.Button(window, text='Batch Treatment', command=lambda: [window.destroy(), batch_treatment()])
    btn_bt.place(x=50, y=200)
    btn_cf = tk.Button(window, text='Continuous Flow Treatment', command=lambda: [window.destroy(), continuous_flow()])
    btn_cf.place(x=230, y=200)

    window.mainloop()


def batch_treatment():
    global days
    global function_choose

    bt = tk.Tk()
    bt.title('Batch Treatment')
    bt.geometry('450x300+400+200')
    function_choose = 2
    days = tk.StringVar()
    days.set(365)

    # user information
    tk.Label(bt, text='Reaction time (day)').place(x=50, y=80)

    entry_day = tk.Entry(bt,textvariable=days)
    entry_day.place(x=240, y=80)

    btn_lan_bt = tk.Button(bt, text='Langmuir isotherms', command=lambda: [bt.destroy(),Langmuir()])
    btn_lan_bt.place(x=50, y=200)

    btn_fre_bt = tk.Button(bt, text='Freundlich isotherms', command=lambda: [bt.destroy(),Freundlich()])
    btn_fre_bt.place(x=200, y=200)

    btn_fre_bt = tk.Button(bt, text='Back', command=lambda: [bt.destroy(),main_window()])
    btn_fre_bt.place(x=360, y=200)

    bt.mainloop()


def continuous_flow():
    global j_in
    global function_choose

    cf = tk.Tk()
    function_choose = 1
    cf.title('Continuous Flow Treatment')
    cf.geometry('450x300+400+200')

    j_in = tk.StringVar()
    j_in.set(0.001)

    # user information
    tk.Label(cf, text='Turnover rate (j) ').place(x=50, y=80)

    entry_j = tk.Entry(cf,textvariable=j_in)
    entry_j.place(x=240, y=80)

    btn_lan_cf = tk.Button(cf, text='Langmuir isotherms', command=lambda: [cf.destroy(),Langmuir()])
    btn_lan_cf.place(x=50, y=200)

    btn_fre_cf = tk.Button(cf, text='Freundlich isotherms', command=lambda: [cf.destroy(),Freundlich()])
    btn_fre_cf.place(x=200, y=200)

    btn_back_cf = tk.Button(cf, text='Back', command=lambda: [cf.destroy(),main_window()])
    btn_back_cf.place(x=360, y=200)

    cf.mainloop()


def Langmuir():
    global Qmax_in
    global KL_in
    Lan = tk.Tk()
    Lan.geometry('450x300+400+200')
    Lan.title('Predicting Removal of Arsenic from Sorbents')
    Qmax_in = tk.StringVar()
    Qmax_in.set(29.44)
    KL_in = tk.StringVar()
    KL_in.set(5.11)

    tk.Label(Lan, text='Maximum adsorption capacity (Qmax) ').place(x=1, y=50)
    tk.Label(Lan, text='Langmuir constant (Kl) ').place(x=1, y=80)

    lan_Qmax = tk.Entry(Lan, textvariable=Qmax_in)
    lan_Qmax.place(x=240, y=50)

    lan_KL = tk.Entry(Lan, textvariable=KL_in)
    lan_KL.place(x=240, y=80)


    btn_Lan = tk.Button(Lan, text='Result', command=lambda: [Lan.destroy(), Result_windows("Langmuir")])
    btn_Lan.place(x=50, y=200)



    btn_Lan_back = tk.Button(Lan, text='Back', command=lambda: [Lan.destroy(), main_window()])
    btn_Lan_back.place(x=230, y=200)


    Lan.mainloop()


def Freundlich():
    global n_in
    global KF_in

    Fre = tk.Tk()
    Fre.geometry('450x300+400+200')
    Fre.title('Predicting Removal of Arsenic from Sorbents')
    KF_in = tk.StringVar()
    KF_in.set(5.10)
    n_in = tk.StringVar()
    n_in.set(2.63)

    tk.Label(Fre, text='constant n ').place(x=1, y=50)
    tk.Label(Fre, text='Freundlich constant (KF) ').place(x=1, y=80)

    Fre_n = tk.Entry(Fre, textvariable=n_in)
    Fre_n.place(x=240, y=50)

    Fre_KF = tk.Entry(Fre, textvariable=KF_in)
    Fre_KF.place(x=240, y=80)

    btn_Lan = tk.Button(Fre, text='Result', command=lambda: [Fre.destroy(), Result_windows()])
    btn_Lan.place(x=60, y=200)
    btn_Lan_back = tk.Button(Fre, text='Back', command=lambda: [Fre.destroy(), main_window()])
    btn_Lan_back.place(x=230, y=200)

    Fre.mainloop()


def plot(time,results_table, C_s_in, C_d_in, isotherms, oder):
    root = tk.Tk()

    if function_choose == 1:
        root.title('As(III) breakthrough curves at influence concentration of %.1d ug/L.' % C_d_in)
    elif function_choose == 2:
        root.title('Batch treatment using the concentration of %.1d ug/L adsorbent for 365 days .' % C_d_in)

    f = Figure(figsize=(8, 5), dpi=100)
    ax = f.add_subplot(111)  # Add sub diagram: 1st in 1 row and 1 column

    # Plotting on the subplot obtained earlier
    ax.plot(time, results_table)
    ax.get_xaxis().get_major_formatter().set_scientific(False)
    ax.set_title('As(III) breakthrough curves with %s adsorption isotherm and %d-order pseudo model' % (isotherms, oder))
    ax.legend(labels=['Cs = %.1d g/L' % C_s_in], loc="center left")
    ax.set_xscale('log')
    ax.set_ylabel('As(aq) (ug/L)')
    ax.set_xlabel('bed volumes treated')

    # Display the drawing to tkinter: create a canvas belonging to root and place the drawing on the canvas
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP,  # Top Alignment
                                fill=tk.BOTH,  # Filling method
                                expand=tk.YES)  # Adjusts with window size

    # The navigation toolbar of matplotlib is shown up (it is not shown by default)
    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=tk.TOP,  # get_tk_widget() gets _tkcanvas
                          fill=tk.BOTH,
                          expand=tk.YES)

    def _quit():
        """This function is called when the exit button is clicked"""
        root.quit()  # Ending the main loop
        root.destroy()  # Destruction window

    # Create a button and bind the above function to it
    button = tk.Button(master=root, text="退出", command=_quit)
    # Button on the bottom
    button.pack(side=tk.BOTTOM)

    # Main loop
    root.mainloop()


def Result_windows(isotherms="Freundlich"):
    C_in = int(Cin.get())
    C_s = int(Cs.get())
    oder = int(var_order.get())
    k_i = float(k_in.get())

    if function_choose == 1:      #CF
        jin = float(j_in.get())
        if isotherms == "Freundlich":
            K_in = float(KF_in.get())
            n_i = float(n_in.get())

            results_table_cf, counter = Result_Processing_cf(j_flow=jin, C_inf=C_in, Cs0=C_s,  ord_in=oder,
                                                                 k_in1=k_i, KF_in1=K_in,n_in1=n_i)

            result = Adsorbents_required(C_in1=C_in, Cs_in=C_s, j_flow=jin, ord_in=oder, isotherms1=isotherms,
                                            k_in1=k_i, KF_in1=K_in, n_in1=n_i)

        elif isotherms == "Langmuir":
            Q_max = float(Qmax_in.get())
            KL = float(KL_in.get())

            results_table_cf, counter = Result_Processing_cf(j_flow=jin, C_inf=C_in, Cs0=C_s,  ord_in=oder,
                                                         isotherms1=isotherms, k_in1=k_i, Qmax_in1=Q_max, KL_in1=KL)

            result = Adsorbents_required(C_in1=C_in, j_flow=jin, Cs_in=C_s,ord_in=oder, isotherms1=isotherms,
                                            k_in1=k_i, Qmax_in1=Q_max, KL_in1=KL)



        result_w = tk.Tk()
        result_w.title('Sorbent efficiency ')
        result_w.geometry('450x300+400+200')
        l = tk.Label(result_w, text='The amount of absorbent required is %.2f kg per year' %result,
                     bg='green',font=('Arial', 12), width=100,
                     height=10)
        l.pack()

        # results_table,C_s_in,C_d_in

        btn_plot_cf = tk.Button(master=result_w, text='plot results', command=lambda: [plot(results_table_cf[0][0],
                                                                                         results_table_cf[0][1], C_s,
                                                                                         C_in, isotherms, oder)])

        btn_plot_cf.place(x=100, y=200)
        btn_Lan_back = tk.Button(result_w, text='Back', command=lambda: [result_w.destroy(), main_window()])
        btn_Lan_back.place(x=230, y=200)
        result_w.mainloop()

    elif function_choose == 2:    # BT
        day = int(days.get())
        result_w = tk.Tk()
        result_w.title('Sorbent efficiency ')
        result_w.geometry('450x300+400+200')
        if isotherms == "Freundlich":
            K_in = float(KF_in.get())
            n_i = float(n_in.get())

            results_Ct, Cs_F, data_collect, results_table_F = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                                         k_i, K_in, n_i)

            result = Min_adsorbent_bt(results_table_F[-1], Cs_F)


        elif isotherms == "Langmuir":
            Q_max = float(Qmax_in.get())
            KL = float(KL_in.get())
            results_Ct, Cs_L, data_collect, results_table_L = Result_Processing_batch(C_in, C_s, day, oder, isotherms, k_i,
                                                         Qmax_in=Q_max, KL_in=KL)

            result = Min_adsorbent_bt(results_table_L[-1], Cs_L)

        l = tk.Label(result_w, text='The amount of absorbent required is %.2f kg per year' % result,
                     bg='green', font=('Arial', 12), width=100,
                     height=10)
        l.pack()
        btn_plot_cf = tk.Button(master=result_w, text='plot results', command=lambda: [plot(data_collect[1:],
                                                                                            results_Ct[4][1:], C_s,
                                                                                            C_in, isotherms, oder)])
        btn_plot_cf.place(x=100, y=200)

        btn_Lan_back = tk.Button(result_w, text='Back', command=lambda: [result_w.destroy(), main_window()])
        btn_Lan_back.place(x=270, y=200)
        result_w.mainloop()


main_window()






