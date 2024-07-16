# Enclosing Elements

No programa JuBEM, Enclosing Elements são utilizados para calcular as integrações singulares fortes através da teoria de deslocamento de corpo rígido em malhas de domínio infinito (em problemas de estática). Para isso, é necessário discretizar um domínio fechado com enclosing elements, que são integrados e cujas entradas são utilizadas no cálculo dos termos singulares fortes.

Após a integração, porém, deve haver uma tratamento para remover os termos de Enclosing Elements do problema, uma vez que não devem ser considerados (se fossem considerados, o domínio não seria infinito). Para isso, foi implementada a função `remove_EE!()`. A seguir esta função será explicada, junto da forma correta de tratar os dados que saem dela.

## Função `remove_EE!()`

