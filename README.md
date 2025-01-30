# Estudo do Método de Diferenças Finitas para Equação da Onda unidimensional

Este repositório explora a resolução numérica da equação da onda unidimensional utilizando o Método de Diferenças Finitas (MDF). O modelo é aplicado em três exemplos distintos, cada um explorando diferentes aspectos da solução numérica, incluindo convergência, simetria e periodicidade. São implementadas diferentes condições de contorno (Dirichlet e Neumann) e analisados os efeitos dos parâmetros numéricos na precisão da solução.

## Descrição do Problema 📝

A equação da onda que estamos resolvendo é dada por:

```math
u_{tt} = c^2u_{xx}, \quad c > 0
```

com condições iniciais:
```math
u(x,0) = Φ(x) \quad \text{e} \quad u_t(x,0) = Ψ(x)
```

e esquema de diferenças finitas:
```math
U^{n+1}_m = c^2λ^2(U^n_{m-1} + U^n_{m+1}) + 2(1-c^2λ^2)U^n_m - U^{n-1}_m
```

## Ferramentas Utilizadas 🔧

Para a implementação numérica, foram utilizadas as seguintes bibliotecas Python:

- **NumPy**: Para operações com arrays e cálculos numéricos
- **Matplotlib.pyplot**: Para visualização dos resultados

## Metodologia 💻

1. **Discretização**: Aplicação do esquema explícito de diferenças finitas
2. **Condições de Contorno**: Implementação de condições de Dirichlet e Neumann
3. **Análise**: Estudo de convergência e erro para diferentes valores de h e λ
