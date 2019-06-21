# Lista5-GuilhermeLyra

A ideia era conseguir resolver o problema do "Shortest Common Superstring", cujo qual é:
```
Dadas N strings, qual é a menor string X que possui as N strings como substrings?

Ex:
Input:  S = {"catg", "ctaagt", "gcta", "ttca", "atgcatc"}
Output: gctaagttcatgcatc
```

E, após isso, comparar o resultado final e o input, utilizando o algoritmo de Sequence Alignment.

Para gerar o input, primeiro foi coletado uma sequência [deste dataset](http://www.vision.ime.usp.br/~jmena/mgwt/datasets/);

Depois, esta sequencia é "estilhaçada" em 12 pedaços de tamanhos aleatórios;

Após isso, tenta-se reconstruir a sequencia original, utilizando o SCS.

# pós-mortem

Embora a ideia fosse ligeiramente boa, as implementações do SCS ou estão:
1. Inconsistente (resultado não ótimo) & Rápido
   
ou

2. Consistente & Lento

# Referencias

> https://fileadmin.cs.lth.se/cs/Personal/Andrzej_Lingas/superstring.pdf

> https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_scs.pdf

> https://www.researchgate.net/profile/Sunil_Kumar_Sahu/post/Does_anyone_have_an_e-book_about_bioinformatics_How_can_i_find_it/attachment/59d6234ac49f478072e995f8/AS%3A272115389927425%401441888772501/download/Bioinformatics+Sequence+and+Genome+Analysis+-+David+W.+Mount.pdf

> https://www8.cs.umu.se/kurser/TDBAfl/VT06/algorithms/BOOK/BOOK5/NODE209.HTM

> https://math.mit.edu/~goemans/18434S06/superstring-lele.pdf

> https://www.geeksforgeeks.org/shortest-superstring-problem-set-2-using-set-cover/

> http://sunny.today/generate-random-integers-with-fixed-sum/

> http://www.vision.ime.usp.br/~jmena/mgwt/datasets/