program ex01_p01
! ---------------------------------------------------------------------------
! Autor: Felipe Baglioli
! Criado em : 21/06/2023
! Disciplina: EAMB-7029 / 2 trim. 2023
! Exercicio 01 - Parte A
! ---------------------------------------------------------------------------
! Avalia tres expressoes com diferentes ordens de operacao, empregando 
! variaveis inteiras e, posteriormente, vairaveis de precisao dupla.
! ---------------------------------------------------------------------------
! VARIAVEIS:
! - a, b, c, d, e, f, g (integer): constantes fornecidas no exercicio para
!   o calculo das expressoes.
!
! - da, db, dc, dd, de, df (real, double precision): constantes acima 
!   convertidas para reais com precisao dupla. A variavel 'g' nao e 
!   convertida, uma vez que se trata de um expoente.
!
! - out1, out2, out3 (integer): respostas das expressoes 1, 2, e 3, 
!   respectivamente, fornecidas no exercicio, com valores inteiros.
!
! - dout1, dout2, dout3 (real, double precision): respostas das expressoes
!   1, 2, e 3, respectivamente, fornecidas no exercicio, com valores reais.
!
! - ndouble (integer): numero kind para precisao dupla.
! ---------------------------------------------------------------------------
! SAIDAS:
! - saida_ex1_a.txt (arquivo): arquivo de saida contendo os resultados das
!   expressoes.
! ---------------------------------------------------------------------------
! Última alteracao em: 23/06/2023

! Declaracao e inicializacao de variaveis:

implicit none

integer :: a = 9, b = 5, c = 8, d = 4, e = 7, f = 3, g = 4
integer :: out1, out2, out3
integer, parameter :: n_double = selected_real_kind(p = 15)

double precision :: da, db, dc, dd, de, df
double precision :: dout1, dout2, dout3

! Abertura do arquivo de saida:

open(unit = 11, file = "saida_ex1_a.txt")

! Calculo das expressoes com valoes inteiros:

out1 = a*b + c*d + e/f**g
out2 = a*(b + c)*d + (e/f)**g
out3 = a*(b + c)*(d + e)/f**g

! Escrita os valores calculados no arquivo de saida:

write(11, *) "Valores Inteiros:"
write(11, *)
write(11, *) "Resultado da expressao 1):", out1
write(11, *) "Resultado da expressao 2):", out2
write(11, *) "Resultado da expressao 3):", out3
write(11, *)

! Conversao dos valores para double precision:

da = real(a, kind = n_double)
db = real(b, kind = n_double)
dc = real(c, kind = n_double)
dd = real(d, kind = n_double)
de = real(e, kind = n_double)
df = real(f, kind = n_double)

! Calculo das expressoes em real double precision:

dout1 = da*db + dc*dd + de/df**g
dout2 = da*(db + dc)*dd + (de/df)**g
dout3 = da*(db + dc)*(dd + de)/df**g

! Escrita os valores calculados no arquivo de saida:

write(11, *) "Valores Reais (precisao dupla):"
write(11, *)
write(11, *) "Resultado da expressao 1):", dout1
write(11, *) "Resultado da expressao 2):", dout2
write(11, *) "Resultado da expressao 3):", dout3
write(11, *)

! Aviso de conclusao do programa:

write(*, *) "Execucao concluida. Por favor, cheque a pasta de saida."


end program ex01_p01