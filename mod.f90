program modch
implicit none
integer(4)op


op=5


if(op==0)then

write(6,*) 'something op=0', mod(op,2)
else if(mod(op,2)==0)then
write(6,*) 'something mod op =0'
else

write(6,*) 'op is odd'
end if





end program modch
