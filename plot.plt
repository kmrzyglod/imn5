# ------------------------------------------------------------------
# Rysowanie gifow
# ------------------------------------------------------------------

reset
set term gif animate size 800,300
set output "anim1.gif"
n=50    #n frames
set pm3d
set view map
set size ratio -1
#set cbrange [0:1]

i=3
load "animacja_zad1.plt"
set output

reset
set term gif animate size 800,300
set output "anim2.gif"
n=50    #n frames
set pm3d
set view map
set size ratio -1
#set cbrange [0:1]

i=3
load "animacja_zad2.plt"
set output

reset
set terminal png size 800,600 enhanced font 'Helvetica,12'
set o 'IntegralZad1.png'
set xl 'czas'
set yl 'wrtosc calki'
#set size ratio -1
set title 'Wartosc calki w kolejnych iteracjach'
#set log x
plot 'zad1_integral.txt' with points 

reset

set terminal png size 800,600 enhanced font 'Helvetica,12'
set o 'IntegralZad2.png'
set xl 'czas'
set yl 'wrtosc calki'
#set size ratio -1
set title 'Wartosc calki w kolejnych iteracjach'
#set log x
plot 'zad2_integral.txt' with points 

reset

set terminal png size 800,600 enhanced font 'Helvetica,12'
set o 'CenterZad1.png'
set xl 'czas'
set yl 'wrtosc calki'
#set size ratio -1
set title 'Wartosc calki w kolejnych iteracjach'
#set log x
plot 'zad1_center.txt' with points 

reset

set terminal png size 800,600 enhanced font 'Helvetica,12'
set o 'CenterZad2.png'
set xl 'czas'
set yl 'wrtosc calki'
#set size ratio -1
set title 'Wartosc calki w kolejnych iteracjach'
#set log x
plot 'zad2_center.txt' with points 



