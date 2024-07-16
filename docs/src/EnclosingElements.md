# Enclosing Elements

No programa JuBEM, Enclosing Elements são utilizados para calcular as integrações singulares fortes através da teoria de deslocamento de corpo rígido em malhas de domínio infinito (em problemas de estática). Para isso, é necessário discretizar um domínio fechado com enclosing elements, que são integrados e cujas entradas são utilizadas no cálculo dos termos singulares fortes.

Após a integração, porém, deve haver uma tratamento para remover os termos de Enclosing Elements do problema, uma vez que não devem ser considerados (se fossem considerados, o domínio não seria infinito). Para isso, foi implementada a função `remove_EE!()`. A seguir esta função será explicada, junto da forma correta de tratar os dados que saem dela.

## Função `remove_EE!()`

Enclosing elements são identificados pelo tipo de condição de contorno `bctype = 0`. O tipo de condição de contorno é armazenado na variável `problem.bctype`. Esta variável possui tamanho nbcx4, onde nbc é o número de condições de contorno diferentes. as colunas contém o número da condição de contorno e o tipo de condição de contorno em cada direção. Uma vez que a condição de enclosing elements não depende da direção, então todas as entradas devem ter valor 0 e qualquer uma pode ser utilizada.

Inicialmente são encontrados os índices das condições de contorno de enclosing elements, chamados `bctypeidxee`. Então são encontrados os índices das tags (atribuídas a cada superfície) que possuem esta condição de contorno, `tagidxee`. Com isso, enfim podem ser encontrados os índices dos elementos que são enclosing elements, `elemidxe`, e os índices dos que não são, `elemidxnonee`.

O que precisamos fazer agora é atualizar as matrizes H e G para conter apenas as entradas de elementos que não são enclosing elements. Para isso, precisamos atualizar as matrizes dos graus de liberdade que são entradas destas matrizes, as matrizes `ID` e `LM`. Além disso, são atualizadas as matrizes de incidência dos elementos físicos e geométricos, `IEN` e `IEN_geo`. 

As matrizes de pontos e nós não são atulizadas, uma vez que a atualização destas iria requerer uma nova ordenação dos nós geométricos na matriz de incidência geométrica `IEN_geo`. Sendo assim, os números globais dos nós são preservados, mas os graus de liberdade dos nós físicos são alterados de forma a remover os graus de liberdade dos enclosing Elements.

Para isso ser possível, a matriz `ID` permanece com o mesmo tamanho que anteriormente, mas as colunas correspondentes aos nós físicos de enclosing elements são zeradas. A matriz `LM` por outro lado, agora tem um tamanho diferente, possuindo número de colunas igual ao número de elementos que não são enclosing elements.

Sendo assim, deve-se levar em consideração que os graus de liberdade devem sempre ser tratados utilizando as matrizes ID e LM.

Por fim, os seguintes valores também são atualizados:


```
    mesh.nelem
    mesh.nnodes
    mesh.ID
    mesh.LM
    mesh.IEN
    mesh.IEN_geo
    mesh.tag
```