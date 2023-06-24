program ex01_p02
! ---------------------------------------------------------------------------
! Autor: Felipe Baglioli
! Criado em : 21/06/2023
! Disciplina: EAMB-7029 / 2 trim. 2023
! Exercicio 01 - Parte B
! ---------------------------------------------------------------------------
! Obtem o resultado de uma expressao extensa com dados fornecidos.
! ---------------------------------------------------------------------------
! VARIAVEIS:
! - a, b, c, d, e, f, g (real): constantes fornecidas no exercicio para
!   o calculo das expressoes.
!
! - y (real): resultado da expressao completa
!
! - p1, p2, p3, p4, p5 (real): etapas intermediarias do calculo de y
!
! ---------------------------------------------------------------------------
! SAIDAS:
! - saida_ex1_b.txt (arquivo): arquivo de saida contendo o resultado final.
! ---------------------------------------------------------------------------
! Última alteracao em: 21/06/2023

! Declaracao e inicializacao de variaveis:

implicit none

real :: a = 9., b = 5., c = 8., d = 4., e = 7., f = 3., g = 4.
real :: y, p1, p2, p3, p4, p5

! Abertura do arquivo de saida:

open(unit = 11, file = "saida_ex1_b.txt")

! Calculo das partes intermediarias:

p1 = 13.1*(a**((2*b)/5.4))

p2 = ((11.4*a*b) + (17*d*c))/(c - (6.7*d))

p3 = g*(f - ((16*d)/(e +(7.8*b*f))))

p4 = ((d*e) + (f**(g/(b**2))))/(b*(c**(e/(g**(1./3.)))))

p5 = ((2*f) + (a**(1./4.)))**(2./3.)

! Resultado final da expressao completa:

y = p1 + p2 + ((p3 + p4)*p5)

! Escrita do resultado no arquivo de saida:

write (11,*) "Resultado final:"
write (11,*)
write (11,*) "y =", y

! Aviso de conclusao do programa:

write(*, *) "Execucao concluida. Por favor, cheque a pasta de saida."

end program ex01_p02