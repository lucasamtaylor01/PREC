import numpy as np
import matplotlib.pyplot as plt


def resolve_eq_onda(a, b, T, h, lambda_val, c=1.0, phi=None, psi=None,
                    cont_esq=None, cont_dir=None,
                    tipo_cont_esq='dirichlet',
                    tipo_cont_dir='dirichlet',
                    ordem_neumann=2):
                      
    # Discretização do domínio 
    M = int((b - a) / h)
    tau = lambda_val * h
    N = int(T / tau)

    x = np.linspace(a, b, M + 1)
    t = np.linspace(0, T, N + 1)
    U = np.zeros((N + 1, M + 1))

    # Condições iniciais
    if phi is not None:
        U[0, :] = phi(x)
    if psi is not None:
        # Primeiro passo temporal (ordem 2)
        U[1, 1:-1] = (c ** 2 * lambda_val ** 2 / 2) * (U[0, :-2] + U[0, 2:]) + \
                     (1 - c ** 2 * lambda_val ** 2) * U[0, 1:-1] + \
                     tau * psi(x[1:-1])
        U[1, 0] = aplica_cond_contorno(U, t, 1, 0, cont_esq, tipo_cont_esq, h, ordem_neumann)
        U[1, -1] = aplica_cond_contorno(U, t, 1, M, cont_dir, tipo_cont_dir, h, ordem_neumann)

    # Esquema explícito para demais passos
    for n in range(1, N):
        U[n + 1, 1:-1] = c ** 2 * lambda_val ** 2 * (U[n, :-2] + U[n, 2:]) + \
                         2 * (1 - c ** 2 * lambda_val ** 2) * U[n, 1:-1] - \
                         U[n - 1, 1:-1]
        U[n + 1, 0] = aplica_cond_contorno(U, t, n + 1, 0, cont_esq, tipo_cont_esq, h, ordem_neumann)
        U[n + 1, -1] = aplica_cond_contorno(U, t, n + 1, M, cont_dir, tipo_cont_dir, h, ordem_neumann)

    return U, x, t


def aplica_cond_contorno(U, t, n, m, func_cont, tipo_cont, h, ordem):
    if tipo_cont.lower() == 'dirichlet':
        return func_cont(t[n])
    elif tipo_cont.lower() == 'neumann':
        if m == 0:  # Contorno esquerdo
            if ordem == 1:
                return U[n, 1] - h * func_cont(t[n])
            return (4 * U[n, 1] - U[n, 2] - 2 * h * func_cont(t[n])) / 3
        else:  # Contorno direito
            if ordem == 1:
                return U[n, -2] + h * func_cont(t[n])
            return (4 * U[n, -2] - U[n, -3] + 2 * h * func_cont(t[n])) / 3


def exemplo1():
    def phi(x): return 2 * np.cos(x)
    def psi(x): return np.zeros_like(x)
    def cont_esq(t): return 0                        # u_x(0,t) = 0
    def cont_dir(t): return 2 * np.cos(1) * np.cos(t)  # u(1,t)
    def sol_exata(x, t): return np.cos(x + t) + np.cos(x - t)

    h_vals = [1/10, 1/20, 1/40]
    ordens = [1, 2]
    resultados = []

    for ordem in ordens:
        print(f"\nOrdem Neumann: {ordem}")
        for h in h_vals:
            U, x, t = resolve_eq_onda(
                a=0, b=1, T=1, h=h, lambda_val=1.0,
                phi=phi, psi=psi,
                cont_esq=cont_esq, cont_dir=cont_dir,
                tipo_cont_esq='neumann', tipo_cont_dir='dirichlet',
                ordem_neumann=ordem
            )

            erro = np.max(np.abs(U[-1, :] - sol_exata(x, t[-1])))
            resultados.append({'h': h, 'ordem': ordem, 'erro': erro})
            print(f"h = {h}, Erro = {erro:.6e}")

            plt.figure(figsize=(10, 6))
            plt.plot(x, U[-1, :], 'b-', label='Numérica')
            plt.plot(x, sol_exata(x, t[-1]), 'r--', label='Exata')
            plt.title(f'Solução em T=1 (h={h}, Ordem={ordem})')
            plt.xlabel('x')
            plt.ylabel('u(x,T)')
            plt.legend()
            plt.grid(True)
            plt.show()

    return resultados


def exemplo2():
    def phi(x): return np.where(np.abs(x) <= 1, 1 - np.abs(x), 0)
    def psi(x): return np.zeros_like(x)
    def cont_esq(t): return 0
    def cont_dir(t): return 0

    h_vals = [1/10, 1/20, 1/40, 1/80]
    lambda_val = 0.95
    resultados = []
    
    for h in h_vals:
        # Desloca condição inicial para x=2
        U, x_calc, _ = resolve_eq_onda(
            a=0, b=4, T=3.8, h=h, lambda_val=lambda_val,
            phi=lambda x: np.where(np.abs(x-2) <= 1, 1 - np.abs(x-2), 0),
            psi=psi,
            cont_esq=cont_esq, cont_dir=cont_dir,
            tipo_cont_esq='dirichlet', tipo_cont_dir='neumann'
        )

        # Verifica simetria em torno de x=2
        idx_2 = np.abs(x_calc - 2).argmin()
        U_esq = U[-1, :idx_2]
        U_dir = U[-1, idx_2:][::-1]
        n = min(len(U_esq), len(U_dir))
        erro_simetria = np.max(np.abs(U_esq[-n:] - U_dir[-n:]))
        resultados.append((h, erro_simetria))

        x_plot = np.linspace(-2, 2, len(x_calc))
        plt.figure(figsize=(10, 6))
        plt.plot(x_plot, U[-1, :])
        plt.title(f'Solução em T=3.8 (h={h})')
        plt.xlabel('x')
        plt.ylabel('u(x,T)')
        plt.grid(True)
        plt.show()

    print("\nAnálise da solução em T=3.8:")
    print("h\t\tErro de simetria")
    print("-" * 30)
    for h, erro in resultados:
        print(f"{h:.6f}\t{erro:.6e}")
        
    return resultados


def exemplo3():
    def phi(x):
        return np.exp(-1000 * (x - 0.5) ** 2) * np.sin(300 * x)
    def psi(x): return np.zeros_like(x)
    def cont(t): return 0

    h = 1/300
    lambda_vals = [1.0, 0.5]
    tempos = [0, 0.25, 2, 10]

    for lambda_val in lambda_vals:
        U, x, t = resolve_eq_onda(
            a=0, b=1, T=10, h=h, lambda_val=lambda_val,
            phi=phi, psi=psi,
            cont_esq=cont, cont_dir=cont,
            tipo_cont_esq='dirichlet', tipo_cont_dir='dirichlet'
        )

        for t_plot in tempos:
            n = int(t_plot / (lambda_val * h))
            plt.figure(figsize=(10, 6))
            plt.plot(x, U[n, :])
            plt.title(f'Solução em t={t_plot} (λ={lambda_val})')
            plt.xlabel('x')
            plt.ylabel('u(x,t)')
            plt.grid(True)
            plt.show()

        # Verifica retorno à condição inicial
        erro_2 = np.max(np.abs(U[int(2 / (lambda_val * h)), :] - U[0, :]))
        erro_10 = np.max(np.abs(U[int(10 / (lambda_val * h)), :] - U[0, :]))
        print(f"\nλ = {lambda_val}")
        print(f"Erro em t=2: {erro_2:.6e}")
        print(f"Erro em t=10: {erro_10:.6e}")


def main():
    while True:
        print("\nEPREC - Equação da Onda\n")
        print("1. Exemplo 1")
        print("2. Exemplo 2")
        print("3. Exemplo 3")
        print("4. Sair")

        try:
            opcao = int(input("\nEscolha um exemplo (1-4): "))
            if opcao == 1:
                exemplo1()
            elif opcao == 2:
                exemplo2()
            elif opcao == 3:
                exemplo3()
            elif opcao == 4:
                print("Encerrando...")
                break
            else:
                print("Opção inválida.")
        except ValueError:
            print("Entrada inválida. Digite um número entre 1 e 4.")


if __name__ == "__main__":
    main()
