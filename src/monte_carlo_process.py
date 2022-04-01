import os
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd


def read_MT_data(path: str):
    output_path = os.path.join(path, 'monte_carlo', 'output')
    cols = ['time', 'T', 'mean_magnetisation_length', 1, 2, 'M']
    df = pd.read_table(output_path, skiprows=8,
                       delim_whitespace=True, names=cols, index_col=False)
    temperatures = df['T']
    magnetizations = df['M']
    return temperatures, magnetizations


def MT_eqation(x: float, Tc: float, beta: float) -> float:
    return (1 - x / Tc) ** beta


def fit_MT_curve(temperatures, magnetizations, threshold: float = .1) -> tuple:
    # pick only temperatures and magnetizations correspond to magnetization higher than treashold (10% of max)
    T_fit = temperatures[magnetizations > threshold].to_numpy()
    M_fit = magnetizations[magnetizations > threshold].to_numpy()
    T_init = T_fit[-1]
    fit_parameters, pcov = curve_fit(
        MT_eqation, T_fit, M_fit, p0=[T_init, 0.4])
    Tc, beta = fit_parameters

    # std of Tc and beta fitting
    Tc_std, beta_std = np.sqrt(np.diag(pcov))
    return Tc, beta, Tc_std, beta_std


def plot_MT_curve(temperature: np.array, magnetization: np.array, Tc: float, beta: float, save_path: str = None) -> None:
    plt.figure(figsize=(8, 5), dpi=150)
    plt.plot(temperature, MT_eqation(temperature, Tc, beta), 'g--',
             label=f'Fit: $(1 - T / {Tc:.0f}) ^ {{{beta:.3f}}}$')
    plt.scatter(temperature, magnetization, label='MC simulation', s=5)
    plt.ylabel('Magnetization')
    plt.xlabel('Temperature, K')
    plt.title('M-T curve', fontweight='bold')
    plt.legend()
    plt.grid(alpha=.4)
    plt.tight_layout()
    if save_path:
        plt.savefig(os.path.join(save_path, 'MT_curve.pdf'))


def form_text_out(Tc: float, beta: float, Tc_std: float, beta_std: float) -> str:
    text = f"""
**********************************************************
*                  Monte-Carlo simulation result         *
**********************************************************

Critical temperature :
    {Tc:.1f} ± {2 * Tc_std:.1f} K

Beta coeff.:
    {beta:.3f} ± {2 * beta_std:.4f}

    """
    return text


def update_output(path: str, text_out: str) -> None:
    out_file_path = os.path.join(path, 'OUTPUT.txt')
    with open(out_file_path, 'a') as f:
        f.write(text_out)
    print(text_out)


def mc_post_process(path: str) -> None:
    temperatures, magnetizations = read_MT_data(path)
    Tc, beta, Tc_std, beta_std = fit_MT_curve(temperatures, magnetizations)
    text_out = form_text_out(Tc, beta, Tc_std, beta_std)
    update_output(path, text_out)
    plot_MT_curve(temperatures, magnetizations, Tc, beta, save_path=path)
    print('Monte-Carlo simulation post-processing complete.')


if __name__ == '__main__':
    mc_post_process(path=os.getcwd())
