import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.special import logsumexp
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import math
import tkinter.messagebox
from matplotlib.figure import Figure
from tkinter import *
from tkinter.ttk import *
import time

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


def Min_adsorbent_bt(results_counter,Cs_col,days,people):
    ads = Cs_col[2:]
    results = [0] * len(ads)
    for h in range(2, len(results_counter)):
        if results_counter[h]==0:
          results_counter[h]=0.0001
        results[h - 2] = (days / results_counter[h]) * ads[h - 2] * 0.001 * (math.ceil(people*7))

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
    global n_in
    global KF_in
    global k_in
    global Qmax_in
    global KL_in

    global n_in_er
    global KF_in_er
    global k_in_er
    global Qmax_in_er
    global KL_in_er

    window = tk.Tk()
    window.title('Predicting Removal of Sorbate from Sorbents')
    window.geometry('600x400+400+200')
    var_order = tk.StringVar()
    var_order.set(1)
    k_in = tk.StringVar()
    k_in.set(0.1111)
    Qmax_in = tk.StringVar()
    Qmax_in.set(29.44)
    KL_in = tk.StringVar()
    KL_in.set(5.11)
    KF_in = tk.StringVar()
    KF_in.set(5.10)
    n_in = tk.StringVar()
    n_in.set(2.63)

    k_in_er = tk.StringVar()
    k_in_er.set(0.015)
    Qmax_in_er = tk.StringVar()
    Qmax_in_er.set(2.99)
    KL_in_er = tk.StringVar()
    KL_in_er.set(0.10)
    KF_in_er = tk.StringVar()
    KF_in_er.set(1.60)
    n_in_er = tk.StringVar()
    n_in_er.set(0.14)

    # user information
    l = tk.Label(window, text='Prediction of the sorbent required to remove sorbate from water',
                 bg='Blue', font=('Arial', 12), width=50,
                 height=2)

    l.pack()  # Label content content area placement, automatic resizing
    tk.Label(window, text='Langmuir Parameters ',
                 bg='Blue', font=('Arial', 10), width=20,
                 height=1).place(x=25, y=70)

    tk.Label(window, text='Maximum adsorption capacity (mg g^-1)  ').place(x=1, y=110)
    tk.Label(window, text='Langmuir constant (Lmg^-1) ').place(x=1, y=140)

    tk.Label(window, text='Freundlich Parameters ',
                 bg='Blue', font=('Arial', 10), width=20,
                 height=1).place(x=25, y=180)

    tk.Label(window, text='constant n ').place(x=1, y=210)
    tk.Label(window, text='Freundlich constant (mg g^-1(mgL^-1)^(-1/n)) ').place(x=1, y=240)

    tk.Label(window, text='Other Parameters ',
             bg='Blue', font=('Arial', 10), width=20,
             height=1).place(x=25, y=270)

    tk.Label(window, text="Rate constant (L g^-1 min^-1)").place(x=1, y=300)
    tk.Label(window, text="Order of pseudo").place(x=1, y=320)

    tk.Label(window, text='Errors',
             bg='Blue', font=('Arial', 10), width=10,
             height=1).place(x=480, y=70)

    lan_Qmax = tk.Entry(window, textvariable=Qmax_in)
    lan_Qmax.place(x=280, y=110)

    lan_Qmax_error = tk.Entry(window, textvariable=Qmax_in_er)
    lan_Qmax_error.place(x=450, y=110)

    lan_KL = tk.Entry(window, textvariable=KL_in)
    lan_KL.place(x=280, y=140)

    lan_KL_error = tk.Entry(window, textvariable=KL_in_er)
    lan_KL_error.place(x=450, y=140)

    Fre_n = tk.Entry(window, textvariable=n_in)
    Fre_n.place(x=280, y=210)

    Fre_n_error = tk.Entry(window, textvariable=n_in_er)
    Fre_n_error.place(x=450, y=210)

    Fre_KF = tk.Entry(window, textvariable=KF_in)
    Fre_KF.place(x=280, y=240)

    Fre_KF_error = tk.Entry(window, textvariable=KF_in_er)
    Fre_KF_error.place(x=450, y=240)




    entry_k = tk.Entry(window, textvariable=k_in)
    entry_k.place(x=280, y=300)

    entry_k_error = tk.Entry(window, textvariable=k_in_er)
    entry_k_error.place(x=450, y=300)

    # Define checkboxes
    entry_O = tk.Checkbutton(window, text="pseudo-second order", variable=var_order, onvalue=2, offvalue=1)
    entry_O.place(x=280, y=320)

    # login and sign up button

    btn_bt = tk.Button(window, text='Langmuir isotherms', command=lambda: [window.destroy(), Transition("Langmuir")])
    btn_bt.place(x=50, y=360)
    btn_cf = tk.Button(window, text='Freundlich isotherms', command=lambda: [window.destroy(), Transition("Freundlich")])
    btn_cf.place(x=340, y=360)

    window.mainloop()


def batch_treatment(isotherms):
    global days
    bt = tk.Tk()
    bt.title('Batch Treatment')
    bt.geometry('450x300+400+200')
    days = tk.StringVar()
    days.set(365)

    # user information
    tk.Label(bt, text='Reaction time (days)').place(x=50, y=80)

    entry_day = tk.Entry(bt,textvariable=days)
    entry_day.place(x=240, y=80)

    btn_re = tk.Button(bt, text='Results', command=lambda: [bt.destroy(), progress(isotherms, 2)])
    btn_re.place(x=100, y=200)

    btn_ba = tk.Button(bt, text='Back', command=lambda: [bt.destroy(), main_window()])
    btn_ba.place(x=270, y=200)

    bt.mainloop()


def continuous_flow(isotherms):
    global j_in

    cf = tk.Tk()

    cf.title('Continuous Flow Treatment')
    cf.geometry('450x300+400+200')

    j_in = tk.StringVar()
    j_in.set(0.001)

    # user information
    tk.Label(cf, text='Turnover rate (min^-1) ').place(x=50, y=80)

    entry_j = tk.Entry(cf, textvariable=j_in)
    entry_j.place(x=240, y=80)

    btn_re = tk.Button(cf, text='Results', command=lambda: [cf.destroy(), progress(isotherms, 1)])
    btn_re.place(x=100, y=200)

    btn_ba = tk.Button(cf, text='Back', command=lambda: [cf.destroy(), main_window()])
    btn_ba.place(x=240, y=200)

    cf.mainloop()


def progress(isotherms, function_choose):

    root = tk.Tk()

    root.title('Loading........')
    root.geometry('450x300+400+200')
    progressbarOne = Progressbar(root, orient=tkinter.HORIZONTAL, length=350, mode='determinate')
    progressbarOne.place(x=40, y=80)

    def preventRepeatedClicks(on):

        if on == 'kai':
            demoBtn = tkinter.Button(root, text="Running", width=20, height=1,
                                     command=lambda: [Result_windows(isotherms,
                                                                     function_choose, progressbarOne, root)

                                                      ])
            demoBtn.place(x=140, y=200)

    preventRepeatedClicks('kai')

    progressbarOne['maximum'] = 100
    progressbarOne['value'] = 0

    root.mainloop()


def Transition(isotherms):
    global Cin
    global Cs
    global people
    Tra = tk.Tk()
    Tra.geometry('450x300+400+200')
    Tra.title('Predicting Removal of Arsenic from Sorbents')

    Cin = tk.StringVar()
    Cin.set(500)
    Cs = tk.StringVar()
    Cs.set(100)
    people = tk.StringVar()
    people.set(5.7)

    tk.Label(Tra, text='Sorbate concentration (ug/L)').place(x=1, y=50)
    tk.Label(Tra, text='Sorbent concentration (g/L)').place(x=1, y=80)
    tk.Label(Tra, text='Average of people per household').place(x=1, y=110)

    entry_sorbent_con = tk.Entry(Tra, textvariable=Cin)
    entry_sorbent_con.place(x=240, y=50)

    entry_sorbate_con = tk.Entry(Tra, textvariable=Cs)
    entry_sorbate_con.place(x=240, y=80)

    entry_sorbate_con = tk.Entry(Tra, textvariable=people)
    entry_sorbate_con.place(x=240, y=110)

    btn_bt = tk.Button(Tra, text='Batch Treatment ', command=lambda: [Tra.destroy(), batch_treatment(isotherms)])
    btn_bt.place(x=50, y=180)

    btn_cf = tk.Button(Tra, text='Continuous Flow Treatment', command=lambda: [Tra.destroy(), continuous_flow(isotherms)])
    btn_cf.place(x=250, y=180)

    btn_back_cf = tk.Button(Tra, text='Back', command=lambda: [Tra.destroy(), main_window()])
    btn_back_cf.place(x=180, y=230)

    Tra.mainloop()


def plot(time, results_table, C_s_in, C_d_in, isotherms, oder,function_choose,num=1):
    root = tk.Tk()

    if function_choose == 1:
        root.title('Breakthrough curves at influence concentration of %.1d ug/L.' % C_d_in)
    elif function_choose == 2:
        root.title('Batch treatment using the concentration of %.1d ug/L adsorbent for 365 days .' % C_d_in)

    if num == 1:
        f = Figure(figsize=(8, 5), dpi=100)
        ax = f.add_subplot(111)  # Add sub diagram: 1st in 1 row and 1 column

        # Plotting on the subplot obtained earlier
        ax.plot(time, results_table)
        ax.get_xaxis().get_major_formatter().set_scientific(False)
        if function_choose == 1:
            ax.set_title(
                'Breakthrough curves with %s adsorption isotherm and %d-order pseudo model' % (isotherms, oder))
            ax.set_xlabel('bed volumes treated')
        elif function_choose == 2:
            ax.set_title(
                'Batch treatment with %s adsorption isotherm and %d-order pseudo model' % (isotherms, oder))
            ax.set_xlabel('time (min)')
        ax.legend(labels=['Cs = %.1d g/L' % C_s_in], loc="center left")
        ax.set_xscale('log')
        ax.set_ylabel('As(aq) (ug/L)')

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
    button = tk.Button(master=root, text="Exit", command=_quit)
    # Button on the bottom
    button.pack(side=tk.BOTTOM)

    # Main loop
    root.mainloop()


def Result_windows(isotherms, function_choose, progressbarOne, root):
    C_in = int(Cin.get())
    C_s = int(Cs.get())
    oder = int(var_order.get())
    k_i = float(k_in.get())

    n_er = float(n_in_er.get())
    KF_er = float(KF_in_er.get())

    k_er = float(k_in_er.get())

    KL_er = float(KL_in_er.get())
    Qmax_er = float(Qmax_in_er.get())

    K_in = float(KF_in.get())
    n_i = float(n_in.get())

    Q_max = float(Qmax_in.get())
    KL = float(KL_in.get())

    peo = float(people.get())

    def but(time, res, iso, fc):
        a = tk.messagebox.askokcancel('Notice',
                                      "Need to plot effect of different reaction times on sorbate concentration "
                                      "in water under current conditions?")

        if a:
            plot(time,
                 res, C_s,
                 C_in, iso, oder, fc)

    if function_choose == 1:      #CF

        jin = float(j_in.get())

        def plot_cf():
            root = tk.Tk()
            root.title('Adsorbents required to treatment')
            root.geometry('850x880+300+0')
            f, ax = plt.subplots(3, 1, figsize=(5, 4), dpi=100)
            k = [(k_i - k_er), k_i, (k_er + k_i)]
            results_table1 = [result_b_k, result, result_u_k]
            # Plotting on the subplot obtained earlier
            ax[0].plot(k, results_table1, color='y', linestyle=':', linewidth=2, marker='o',
                       markerfacecolor='blue', markersize=8)

            for a, b in zip(k, results_table1):
                ax[0].text(round(a, 3), round(b, 3), round(b, 3))

            ax[0].set_title('Effect of rate constants on adsorbent')
            ax[0].set_xlabel('rate constant ')
            ax[0].set_ylabel('sorbent need ')

            if isotherms == "Freundlich":
                KF1 = [ (K_in - KF_er), K_in, (K_in+KF_er)]
                results_table2 = [result_b_KF, result, result_u_KF]
                # Plotting on the subplot obtained earlier
                ax[1].plot(KF1, results_table2, color='y', linestyle=':', linewidth=2, marker='o',
                           markerfacecolor='blue', markersize=8)

                for a, b in zip(KF1, results_table2):
                    ax[1].text(round(a, 3), round(b, 3), round(b, 3))

                ax[1].set_title('Effect of Freundlich parameters on adsorbent')
                ax[1].set_xlabel('Freundlich parameters ')
                ax[1].set_ylabel('sorbent need ')

                n = [(n_i - n_er), n_i, (n_i + n_er)]
                results_table3 = [result_b_n, result, result_u_n]
                # Plotting on the subplot obtained earlier
                ax[2].plot(n, results_table3, color='y', linestyle=':', linewidth=2, marker='o',
                           markerfacecolor='blue', markersize=8)

                for a, b in zip(n, results_table3):
                    ax[2].text(round(a, 3), round(b, 3), round(b, 3))

                ax[2].set_title('Effect of parameters n on adsorbent')
                ax[2].set_xlabel('parameters n ')
                ax[2].set_ylabel('sorbent need ')

            elif isotherms == "Langmuir":
                KL_X = [ (KL - KL_er), KL, (KL+KL_er)]
                results_table2 = [result_b_KL, result, result_u_KL]
                # Plotting on the subplot obtained earlier
                ax[1].plot(KL_X, results_table2, color='y', linestyle=':', linewidth=2, marker='o',
                           markerfacecolor='blue', markersize=8)

                for a, b in zip(KL_X, results_table2):
                    ax[1].text(round(a, 3), round(b, 3), round(b, 3))
                ax[1].set_title('Effect of Langmuir parameters on adsorbent')
                ax[1].set_xlabel('Langmuir parameters')
                ax[1].set_ylabel('sorbent need ')
                Qmax1 = [(Q_max - Qmax_er), Q_max, (Q_max + Qmax_er)]
                results_table3 = [result_b_Q, result, result_u_Q]
                # Plotting on the subplot obtained earlier
                ax[2].plot(Qmax1, results_table3, color='y', linestyle=':', linewidth=2, marker='o',
                           markerfacecolor='blue', markersize=8)

                for a, b in zip(Qmax1, results_table3):
                    ax[2].text(round(a, 3), round(b, 3), round(b, 3))

                ax[2].set_title('Effect of Maximum adsorption capacity on adsorbent')
                ax[2].set_xlabel('Maximum adsorption capacity')
                ax[2].set_ylabel('sorbent need ')

            plt.tight_layout()
            # Display the drawing to tkinter: create a canvas belonging to root and place the drawing on the canvas
            canvas_spice = FigureCanvasTkAgg(f, root)
            canvas_spice.draw()
            canvas_spice.get_tk_widget().pack(side=tkinter.TOP,
                                              fill=tkinter.BOTH,
                                              expand=tkinter.YES)

            # The navigation toolbar of matplotlib is shown up (it is not shown by default)
            toolbar = NavigationToolbar2Tk(canvas_spice, root)
            toolbar.update()
            canvas_spice._tkcanvas.pack(side=tkinter.TOP,  # get_tk_widget()得到的就是_tkcanvas
                                        fill=tkinter.BOTH,
                                        expand=tkinter.YES)

            def _quit():
                """This function is called when the exit button is clicked"""
                root.quit()  # Ending the main loop
                root.destroy()  # Destruction window

                # Create a button and bind the above function to it

            button = tk.Button(master=root, text="Exit", command=_quit)
            # Button on the bottom
            button.pack(side=tkinter.BOTTOM)

            # Main loop
            root.mainloop()

        if isotherms == "Freundlich":
            results_table_cf, counter = Result_Processing_cf(j_flow=jin, C_inf=C_in, Cs0=C_s,  ord_in=oder,
                                                             k_in1=k_i, KF_in1=K_in, n_in1=n_i)

            progressbarOne['value'] +=12.5

            root.update()

            result = Adsorbents_required(C_in1=C_in, Cs_in=C_s, j_flow=jin, ord_in=oder, isotherms1=isotherms,
                                         k_in1=k_i, KF_in1=K_in, n_in1=n_i)

            progressbarOne['value'] += 12.5

            root.update()

            result_u_KF = Adsorbents_required(C_in1=C_in, Cs_in=C_s, j_flow=jin, ord_in=oder, isotherms1=isotherms,
                                              k_in1=k_i, KF_in1=(K_in+KF_er), n_in1=n_i)

            progressbarOne['value'] += 12.5

            root.update()

            result_b_KF = Adsorbents_required(C_in1=C_in, Cs_in=C_s, j_flow=jin, ord_in=oder, isotherms1=isotherms,
                                              k_in1=k_i, KF_in1=(K_in-KF_er), n_in1=n_i)

            progressbarOne['value'] += 12.5

            root.update()

            result_u_k = Adsorbents_required(C_in1=C_in, Cs_in=C_s,
                                             j_flow=jin, ord_in=oder, isotherms1=isotherms,
                                             k_in1=(k_i+k_er), KF_in1=K_in, n_in1=n_i)

            progressbarOne['value'] += 12.5

            root.update()

            result_b_k = Adsorbents_required(C_in1=C_in, Cs_in=C_s, j_flow=jin, ord_in=oder, isotherms1=isotherms,
                                             k_in1=(k_i-k_er), KF_in1=K_in, n_in1=n_i)

            progressbarOne['value'] += 12.5

            root.update()

            result_u_n = Adsorbents_required(C_in1=C_in, Cs_in=C_s, j_flow=jin, ord_in=oder, isotherms1=isotherms,
                                             k_in1=k_i, KF_in1=K_in, n_in1=n_i+n_er)

            progressbarOne['value'] += 12.5

            root.update()

            result_b_n = Adsorbents_required(C_in1=C_in, Cs_in=C_s, j_flow=jin, ord_in=oder, isotherms1=isotherms,
                                             k_in1=k_i, KF_in1=K_in, n_in1=n_i-n_er)

            progressbarOne['value'] += 12.5

            root.update()

            time.sleep(1.5)
            root.destroy()
            result_w = tk.Tk()
            result_w.title('Sorbent efficiency ')
            result_w.geometry('550x300+400+200')

            tk.Label(result_w, text='Effect of rate constants on adsorbent ',
                     bg='green', font=('Arial', 10), width=50,
                     height=2).place(x=1, y=10)

            tk.Label(result_w,
                     text='The amount of absorbent required per year ranges from %.2f kg to %.2f kg' % (
                     result_u_k, result_b_k)
                     ).place(x=1, y=50)

            tk.Label(result_w, text='Effect of parameter n on adsorbent ',
                     bg='green', font=('Arial', 10), width=50,
                     height=2).place(x=1, y=80)

            tk.Label(result_w,
                     text='The amount of absorbent required per year ranges from %.2f kg to %.2f kg' % (
                         result_u_n, result_b_n)
                     ).place(x=1, y=120)

            tk.Label(result_w, text='Effect of Freundlich parameters on adsorbent ',
                     bg='green', font=('Arial', 10), width=50,
                     height=2).place(x=1, y=150)

            tk.Label(result_w,
                     text='The amount of absorbent required per year ranges from %.2f kg to %.2f kg' % (
                         result_u_KF, result_b_KF)
                     ).place(x=1, y=190)

            btn_plot1_cf = tk.Button(master=result_w, text='plot (time) ', command=lambda: [
                                                                                            but(results_table_cf[0][0],
                                                                                                results_table_cf[0][1],
                                                                                                isotherms, 1)])

            btn_plot2_cf = tk.Button(master=result_w, text='plot results', command=plot_cf)

            btn_plot1_cf.place(x=230, y=250)
            btn_plot2_cf.place(x=100, y=250)
            btn_Lan_back = tk.Button(result_w, text='Back',
                                     command=lambda: [ result_w.destroy(), main_window()])
            btn_Lan_back.place(x=350, y=250)

            result_w.mainloop()


        elif isotherms == "Langmuir":


            results_table_cf, counter = Result_Processing_cf(j_flow=jin, C_inf=C_in, Cs0=C_s,  ord_in=oder,
                                                             isotherms1=isotherms, k_in1=k_i, Qmax_in1=Q_max, KL_in1=KL)

            progressbarOne['value'] += 12.5

            root.update()

            result = Adsorbents_required(C_in1=C_in, j_flow=jin, Cs_in=C_s, ord_in=oder, isotherms1=isotherms,
                                         k_in1=k_i, Qmax_in1=Q_max, KL_in1=KL)

            progressbarOne['value'] += 12.5

            root.update()

            result_u_k = Adsorbents_required(C_in1=C_in, j_flow=jin, Cs_in=C_s, ord_in=oder, isotherms1=isotherms,
                                             k_in1=(k_i+k_er), Qmax_in1=Q_max, KL_in1=KL)

            progressbarOne['value'] += 12.5

            root.update()

            result_b_k = Adsorbents_required(C_in1=C_in, j_flow=jin, Cs_in=C_s, ord_in=oder, isotherms1=isotherms,
                                             k_in1=(k_i-k_er), Qmax_in1=Q_max, KL_in1=KL)

            progressbarOne['value'] += 12.5

            root.update()

            result_u_KL = Adsorbents_required(C_in1=C_in, j_flow=jin, Cs_in=C_s, ord_in=oder, isotherms1=isotherms,
                                              k_in1=k_i, Qmax_in1=Q_max, KL_in1=KL+KL_er)

            progressbarOne['value'] += 12.5

            root.update()

            result_b_KL = Adsorbents_required(C_in1=C_in, j_flow=jin, Cs_in=C_s, ord_in=oder, isotherms1=isotherms,
                                              k_in1=k_i, Qmax_in1=Q_max, KL_in1=KL-KL_er)

            progressbarOne['value'] += 12.5

            root.update()

            result_u_Q = Adsorbents_required(C_in1=C_in, j_flow=jin, Cs_in=C_s, ord_in=oder, isotherms1=isotherms,
                                             k_in1=k_i, Qmax_in1=Q_max+Qmax_er, KL_in1=KL)

            progressbarOne['value'] += 12.5

            root.update()

            result_b_Q = Adsorbents_required(C_in1=C_in, j_flow=jin, Cs_in=C_s, ord_in=oder, isotherms1=isotherms,
                                             k_in1=k_i, Qmax_in1=Q_max-Qmax_er, KL_in1=KL)

            progressbarOne['value'] += 12.5

            root.update()

            time.sleep(1.5)
            root.destroy()
            result_w = tk.Tk()
            result_w.title('Sorbent efficiency ')
            result_w.geometry('550x300+400+200')

            tk.Label(result_w, text='Effect of rate constants on adsorbent ',
                     bg='green', font=('Arial', 10), width=50,
                     height=2).place(x=1, y=10)

            tk.Label(result_w,
                     text='The amount of absorbent required per year ranges from %.2f kg to %.2f kg' %(result_u_k, result_b_k)
                     ).place(x=1, y=50)

            tk.Label(result_w, text='Effect of Langmuir parameters on adsorbent ',
                     bg='green', font=('Arial', 10), width=50,
                     height=2).place(x=1, y=80)

            tk.Label(result_w,
                     text='The amount of absorbent required per year ranges from %.2f kg to %.2f kg' % (
                      result_u_KL, result_b_KL)
                     ).place(x=1, y=120)

            tk.Label(result_w, text='Effect of Maximum adsorption capacity on adsorbent ',
                     bg='green', font=('Arial', 10), width=50,
                     height=2).place(x=1, y=150)

            tk.Label(result_w,
                     text='The amount of absorbent required per year ranges from %.2f kg to %.2f kg' % (
                      result_u_Q, result_b_Q)
                     ).place(x=1, y=190)

            btn_plot1_cf = tk.Button(master=result_w, text='plot (time) ', command=lambda: [
                                                                                            but(results_table_cf[0][0],
                                                                                                results_table_cf[0][1],
                                                                                                isotherms, 1)])

            btn_plot2_cf = tk.Button(master=result_w, text='plot results', command=lambda: plot_cf)

            btn_plot1_cf.place(x=230, y=250)
            btn_plot2_cf.place(x=100, y=250)
            btn_Lan_back = tk.Button(result_w, text='Back',
                                     command=lambda: [result_w.destroy(), main_window()])
            btn_Lan_back.place(x=350, y=250)

            result_w.mainloop()

    elif function_choose == 2:    # BT
        day = int(days.get())
        if isotherms == "Freundlich":

            # k_in=0.1111, KF_in=5.10, n_in=2.63
            results_Ct, Cs_F, data_collect, results_table_F = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                                                                      k_i, K_in, n_i)

            progressbarOne['value'] += 7.143

            root.update()

            results_Ct_n_u, Cs_F_n_u, data_collect_n_u, results_table_F_n_u \
                = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                          k_i,
                                          n_in=n_i + n_er, KF_in=K_in)

            progressbarOne['value'] += 7.143

            root.update()

            results_Ct_n_b, Cs_F_n_b, data_collect_n_b, results_table_F_n_b \
                = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                          k_i,
                                          n_in=n_i - n_er, KF_in=K_in)

            progressbarOne['value'] += 7.143

            root.update()

            results_Ct_k_u, Cs_F_k_u, data_collect_k_u, results_table_F_k_u \
                = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                          k_i + k_er,
                                          n_in=n_i, KF_in=K_in)

            progressbarOne['value'] += 7.143

            root.update()

            results_Ct_k_b, Cs_F_k_b, data_collect_k_b, results_table_F_k_b \
                = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                          k_i - k_er,
                                          n_in=n_i, KF_in=K_in)

            progressbarOne['value'] += 7.143

            root.update()

            results_Ct_KF_u, Cs_F_KF_u, data_collect_KF_u, results_table_F_KF_u \
                = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                          k_i,
                                          n_in=n_i, KF_in=K_in+KF_er)

            progressbarOne['value'] += 7.143

            root.update()

            results_Ct_KF_b, Cs_F_KF_b, data_collect_KF_b, results_table_F_KF_b \
                = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                          k_i,
                                          n_in=n_i, KL_in=K_in - KL_er)

            progressbarOne['value'] += 7.143

            root.update()

            result = Min_adsorbent_bt(results_table_F[-1], Cs_F, day, peo)

            progressbarOne['value'] += 7.143

            root.update()

            result_n_u = Min_adsorbent_bt(results_table_F_n_u[-1], Cs_F_n_u, day, peo)

            progressbarOne['value'] += 7.143

            root.update()

            result_n_b = Min_adsorbent_bt(results_table_F_n_b[-1], Cs_F_n_b, day, peo)

            progressbarOne['value'] += 7.143

            root.update()

            result_k_u = Min_adsorbent_bt(results_table_F_k_u[-1], Cs_F_k_u, day, peo)

            progressbarOne['value'] += 7.143

            root.update()

            result_k_b = Min_adsorbent_bt(results_table_F_k_b[-1], Cs_F_k_b, day, peo)

            progressbarOne['value'] += 7.143

            root.update()

            result_KF_u = Min_adsorbent_bt(results_table_F_KF_u[-1], Cs_F_KF_u, day, peo)

            progressbarOne['value'] += 7.143

            root.update()

            result_KF_b = Min_adsorbent_bt(results_table_F_KF_b[-1], Cs_F_KF_b, day, peo)

            progressbarOne['value'] += 7.143

            root.update()

            time.sleep(0.5)
            root.destroy()
            result_w = tk.Tk()
            result_w.title('Sorbent efficiency ')
            result_w.geometry('550x330+400+200')

            tk.Label(result_w, text='Effect of rate constants on adsorbent ',
                     bg='green', font=('Arial', 10), width=50,
                     height=2).place(x=1, y=10)

            tk.Label(result_w,
                     text='The amount of absorbent required per %d days ranges from %.2f kg to %.2f kg' % (day,
                                                                                                           result_k_u,
                                                                                                           result_k_b)
                     ).place(x=1, y=50)

            tk.Label(result_w, text='Effect of Freundlich parameters on adsorbent ',
                     bg='green', font=('Arial', 10), width=50,
                     height=2).place(x=1, y=80)

            tk.Label(result_w,
                     text='The amount of absorbent required per %d days ranges from %.2f kg to %.2f kg' % (day,
                                                                                                           result_KF_u,
                                                                                                           result_KF_b)
                     ).place(x=1, y=120)

            tk.Label(result_w, text='Effect of parameters n on adsorbent ',
                     bg='green', font=('Arial', 10), width=50,
                     height=2).place(x=1, y=150)

            tk.Label(result_w,
                     text='The amount of absorbent required per %d days ranges from %.2f kg to %.2f kg' % (day,
                                                                                                           result_n_u,
                                                                                                           result_n_b)
                     ).place(x=1, y=190)

        elif isotherms == "Langmuir":

            results_Ct, Cs_L, data_collect, results_table_L = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                                                                      k_i,
                                                                                      Qmax_in=Q_max, KL_in=KL)
            progressbarOne['value'] += 7.143

            root.update()

            results_Ct_Q_u, Cs_L_Q_u, data_collect_Q_u, results_table_L_Q_u \
                = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                          k_i,
                                          Qmax_in=Q_max+Qmax_er, KL_in=KL)

            progressbarOne['value'] += 7.143

            root.update()

            results_Ct_Q_b, Cs_L_Q_b, data_collect_Q_b, results_table_L_Q_b \
                = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                          k_i,
                                          Qmax_in=Q_max-Qmax_er, KL_in=KL)

            progressbarOne['value'] += 7.143

            root.update()

            results_Ct_k_u, Cs_L_k_u, data_collect_k_u, results_table_L_k_u\
                = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                          k_i+k_er,
                                          Qmax_in=Q_max, KL_in=KL)

            progressbarOne['value'] += 7.143

            root.update()

            results_Ct_k_b, Cs_L_k_b, data_collect_k_b, results_table_L_k_b \
                = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                          k_i-k_er,
                                          Qmax_in=Q_max, KL_in=KL)

            progressbarOne['value'] += 7.143

            root.update()

            results_Ct_KL_u, Cs_L_KL_u, data_collect_KL_u, results_table_L_KL_u \
                = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                          k_i,
                                          Qmax_in=Q_max, KL_in=KL)

            progressbarOne['value'] += 7.143

            root.update()

            results_Ct_KL_b, Cs_L_KL_b, data_collect_KL_b, results_table_L_KL_b \
                = Result_Processing_batch(C_in, C_s, day, oder, isotherms,
                                          k_i,
                                          Qmax_in=Q_max, KL_in=KL-KL_er)

            progressbarOne['value'] += 7.143

            root.update()

            result = Min_adsorbent_bt(results_table_L[-1], Cs_L, day, peo)

            progressbarOne['value'] += 7.143

            root.update()

            result_Q_u = Min_adsorbent_bt(results_table_L_Q_u[-1], Cs_L_Q_u, day, peo)

            progressbarOne['value'] += 7.143

            root.update()

            result_Q_b = Min_adsorbent_bt(results_table_L_Q_b[-1], Cs_L_Q_b, day, peo)

            progressbarOne['value'] += 7.143

            root.update()


            result_k_u = Min_adsorbent_bt(results_table_L_k_u[-1], Cs_L_k_u, day, peo)

            progressbarOne['value'] += 7.143

            root.update()

            result_k_b = Min_adsorbent_bt(results_table_L_k_b[-1], Cs_L_k_b, day, peo)

            progressbarOne['value'] += 7.143

            root.update()

            result_KL_u = Min_adsorbent_bt(results_table_L_KL_u[-1], Cs_L_KL_u, day, peo)

            progressbarOne['value'] += 7.143

            root.update()

            result_KL_b = Min_adsorbent_bt(results_table_L_KL_b[-1], Cs_L_KL_b, day, peo)

            progressbarOne['value'] += 7.143

            root.update()

            time.sleep(0.5)
            root.destroy()
            result_w = tk.Tk()
            result_w.title('Sorbent efficiency ')
            result_w.geometry('550x330+400+200')

            tk.Label(result_w, text='Effect of rate constants on adsorbent ',
                     bg='green', font=('Arial', 10), width=50,
                     height=2).place(x=1, y=10)

            tk.Label(result_w,
                     text='The amount of absorbent required per %d days ranges from %.2f kg to %.2f kg' % (day,
                                                                                                           result_k_u,
                                                                                                           result_k_b)
                     ).place(x=1, y=50)

            tk.Label(result_w, text='Effect of Langmuir parameters on adsorbent ',
                     bg='green', font=('Arial', 10), width=50,
                     height=2).place(x=1, y=80)

            tk.Label(result_w,
                     text='The amount of absorbent required per %d days ranges from %.2f kg to %.2f kg' % (day,
                                                                                                           result_KL_u,
                                                                                                           result_KL_b)
                     ).place(x=1, y=120)

            tk.Label(result_w, text='Effect of Maximum adsorption capacity on adsorbent ',
                     bg='green', font=('Arial', 10), width=50,
                     height=2).place(x=1, y=150)

            tk.Label(result_w,
                     text='The amount of absorbent required per %d days ranges from %.2f kg to %.2f kg' % (day,
                                                                                                           result_Q_u,
                                                                                                           result_Q_b)
                     ).place(x=1, y=190)

        def plot_bt():
            root = tk.Tk()
            root.title('Adsorbents required to treatment')
            root.geometry('850x880+300+0')
            f, ax = plt.subplots(3, 1, figsize=(5, 4), dpi=100)
            k = [(k_i - k_er), k_i, (k_er + k_i)]
            results_table1 = [result_k_b, result, result_k_u]
            # Plotting on the subplot obtained earlier
            ax[0].plot(k, results_table1, color='y', linestyle=':', linewidth=2, marker='o',
                       markerfacecolor='blue', markersize=8)

            for a, b in zip(k, results_table1):
                ax[0].text(round(a, 3), round(b, 3), round(b, 3))

            ax[0].set_title('Effect of rate constants on adsorbent')
            ax[0].set_xlabel('rate constant ')
            ax[0].set_ylabel('sorbent need ')

            if isotherms == "Freundlich":
                KF1 = [(K_in - KF_er), K_in, (K_in + KF_er)]
                results_table2 = [result_KF_b, result, result_KF_u]
                # Plotting on the subplot obtained earlier
                ax[1].plot(KF1, results_table2, color='y', linestyle=':', linewidth=2, marker='o',
                           markerfacecolor='blue', markersize=8)

                for a, b in zip(KF1, results_table2):
                    ax[1].text(round(a, 3), round(b, 3), round(b, 3))

                ax[1].set_title('Effect of Freundlich parameters on adsorbent')
                ax[1].set_xlabel('Freundlich parameters ')
                ax[1].set_ylabel('sorbent need ')

                n = [(n_i - n_er), n_i, (n_i + n_er)]
                results_table3 = [result_n_b, result, result_n_u]
                # Plotting on the subplot obtained earlier
                ax[2].plot(n, results_table3, color='y', linestyle=':', linewidth=2, marker='o',
                           markerfacecolor='blue', markersize=8)

                for a, b in zip(n, results_table3):
                    ax[2].text(round(a, 3), round(b, 3), round(b, 3))

                ax[2].set_title('Effect of parameters n on adsorbent')

                ax[2].set_xlabel('parameters n ')
                ax[2].set_ylabel('sorbent need ')

            elif isotherms == "Langmuir":
                KL_X = [(KL - KL_er), KL, (KL + KL_er)]
                results_table2 = [result_KL_b, result, result_KL_u]
                # Plotting on the subplot obtained earlier
                ax[1].plot(KL_X, results_table2, color='y', linestyle=':', linewidth=2, marker='o',
                           markerfacecolor='blue', markersize=8)

                for a, b in zip(KL_X, results_table2):
                    ax[1].text(round(a, 3), round(b, 3), round(b, 3))
                ax[1].set_title('Effect of Langmuir parameters on adsorbent')
                ax[1].set_xlabel('Langmuir parameters')
                ax[1].set_ylabel('sorbent need ')
                Qmax1 = [(Q_max - Qmax_er), Q_max, (Q_max + Qmax_er)]
                results_table3 = [result_Q_b, result, result_Q_u]
                # Plotting on the subplot obtained earlier
                ax[2].plot(Qmax1, results_table3, color='y', linestyle=':', linewidth=2, marker='o',
                           markerfacecolor='blue', markersize=8)

                for a, b in zip(Qmax1, results_table3):
                    ax[2].text(round(a, 3), round(b, 3), round(b, 3))

                ax[2].set_title('Effect of Maximum adsorption capacity on adsorbent')
                ax[2].set_xlabel('Maximum adsorption capacity')
                ax[2].set_ylabel('sorbent need ')

            plt.tight_layout()
            # Display the drawing to tkinter: create a canvas belonging to root and place the drawing on the canvas
            canvas_spice = FigureCanvasTkAgg(f, root)
            canvas_spice.draw()
            canvas_spice.get_tk_widget().pack(side=tkinter.TOP,
                                              fill=tkinter.BOTH,
                                              expand=tkinter.YES)

            # The navigation toolbar of matplotlib is shown up (it is not shown by default)
            toolbar = NavigationToolbar2Tk(canvas_spice, root)
            toolbar.update()
            canvas_spice._tkcanvas.pack(side=tkinter.TOP,  # get_tk_widget()得到的就是_tkcanvas
                                        fill=tkinter.BOTH,
                                        expand=tkinter.YES)

            def _quit():
                """This function is called when the exit button is clicked"""
                root.quit()  # Ending the main loop
                root.destroy()  # Destruction window

                # Create a button and bind the above function to it

            button = tk.Button(master=root, text="Exit", command=_quit)
            # Button on the bottom
            button.pack(side=tkinter.BOTTOM)

            # Main loop
            root.mainloop()

        btn_plot1_cf = tk.Button(master=result_w, text='plot(time)',
                                 command=lambda: [but(data_collect[1:],
                                                  results_Ct[4][1:],
                                                  isotherms, 2)])

        btn_plot2_cf = tk.Button(master=result_w, text='plot results', command=plot_bt)

        btn_plot1_cf.place(x=230, y=250)
        btn_plot2_cf.place(x=100, y=250)

        btn_Lan_back = tk.Button(result_w, text='Back', command=lambda: [result_w.destroy(), main_window()])
        btn_Lan_back.place(x=350, y=250)
        result_w.mainloop()


main_window()






