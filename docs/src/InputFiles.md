# Esboço da documentação do novo código de JuBEM

## Arquivos de entrada

Na nova versão do JuBEM a intenção é separar os arquivos de entrada em dois, um que contém apenas os dados de malha, sendo o arquivo que vem direto do Gmsh em formato `.msh`, e outro que contém todas as outras informações necessárias para definição do problema e dos parâmetros de solução em formato `.prob`.

Arquivos de malha gerados pelo gmsh, no momento, devem diferenciar as superfícies do problema ao atribuir superfícies físicas a estes. Esta diferenciação é utilizada dentro do programa para definir as condições de contorno das diferentes superfícies. Existe uma intenção em colocar outros tipos de entidades físicas também, como nós, arestas e etc, mas no momento isto ainda não é feito. A diferenciação é realizada por superfícies.

## Arquivo `.prob`

Arquivos `.prob` podem ter os seguintes blocos:

- **ProblemName** (Obrigatório): Contém o nome do problema
- **Material** (Obrigatório): Contém as propriedades de cada material do problema
- **MeshType** (Obrigatório): Contém informações da ordem do elemento e do valor da descontinuidade.
- **Frequencies**: Contém as faixas de frequência de solução do problema dinâmico. (obrigatório para problema dinâmico)
- **BoundaryConditions** (Obrigatório): Contém as informações das diferentes condições de contorno, sendo elas o tipo da condição de contorno e seu valor.
- **TagInformation** (Obrigatório): Contém a associação entre o número da tag física proveniente do gmsh com os valores de condições de contorno e número do material.
- **Forces**: Contém as informações das forças atuando sobre o corpo rígido. (Obrigatório para corpo rígido)

### Descrições dos blocos dos arquivo `.prob`

#### Bloco `TagInformation`

O bloco TagInformation contém as seguintes informações:

```
Ns
idx_s1 bcidx_s1 matidx_s1 bodyidx_s1
idx_s2 bcidx_s2 matidx_s2 bodyidx_s2
...
idx_sNs bcidx_sNs matidx_sNs bodyidx_sNs
```

`Ns` é o número de superfícies físicas, `idx_si` é o índice da superfície física, `bcidx_si` é o índice da condição de contorno da superfície física i, `matidx_si` é o índice do material da superfície física i, `bodyidx_si` é o índice do corpo da superfície física i.

#### Bloco `BoundaryConditions`

O bloco BoundaryConditions contém as seguintes informações:

```
Nbc
idx_bc1 tipo_bcx valor_bcx tipo_bcy valor_bcy tipo_bcz valor_bcz
idx_bc2 tipo_bcx valor_bcx tipo_bcy valor_bcy tipo_bcz valor_bcz
...
idx_bcNs tipo_bcx valor_bcx tipo_bcy valor_bcy tipo_bcz valor_bcz
```

`Nbc` é o número de condições de contorno, `idx_bci` é o índice da condição de contorno e em seguida é passada uma string que contém o tipo de condição de contorno em cada direção e seu valor. Condições de contorno suportadas são:

- `u`: Deslocamento (Dirichlet)
- `t`: Tração superficial (Neumann)
- `rb`: Corpo rígido (Neste caso os valores são a posição do centróide)
- `ee`: Enclosing Element (valores não tem importancia, mas devem ser passados, podem ser 0)
- `i`: Interface (valores devem ser 1 e -1 para cada lado da interface)

No programa JuBEM, cada uma destas será atribuída um valor numérico na variável `problem.bctype`, a relação dos valores é a seguinte:

- `u`: 1
- `t`: 2
- `rb`: 2+(idxrb)
- `ee`: 0
- `i`: -1