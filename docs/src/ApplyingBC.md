# Aplicação das condições de contorno

## Aplicação das condições de contorno de solo homogêneo com N corpos rígidos (função `apply_BC`)

No momento esta função funciona apenas para solos homogêneos. A intenção é que no futuro sirva também para heterogêneos

Esboço:

- Definição do número de nós por elemento `nnel`
- Criação de uma nova matriz LM
    - Esta matriz tem tamanho nnel,nel,3
    - Contém o grau de liberdade do nó local x elemento x direção
- Definição de nu (número de condições de contorno de u)
- Definição de nt (número de condições de contorno de t)
- Definição dos índices dos elementos que possuem condições de contorno u e t (não vai funcionar diretamente com a nova organização dos dados)
- Encontra os índices dos elementos que são flexíveis
- Conta quantos corpos rígidos existem `nrrb`
- Cria um vetor que armazena quantos elementos existem em cada corpo rígido.
- Cria vetor que contém os índices dos elementos rígidos
- Cria vetores abstratos que vão conter os índices dos elementos de cada CR, e as matrizes C e D de cada CR
- 