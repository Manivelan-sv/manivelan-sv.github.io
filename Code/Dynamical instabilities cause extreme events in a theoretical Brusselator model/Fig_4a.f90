! Dynamical instabilities cause extreme events in a theoretical Brusselator model.
! Manivelan, S. V., Sabarathinam, S., Thamilmaran, K., & Manimehan, I.
! Chaos, Solitons & Fractals, 180, 114582. (https://doi.org/10.1016/j.chaos.2024.114582)

implicit real*8(a-h,o-z)
dimension xxl(10000000), xl(10000000), ey(4000000)
Character*20 in, out

! Prompt user to enter input (Data) filename
Print *, 'Enter input (Data) filename :'
Read(*, 11) in

! Prompt user to enter output (PDF) filename
Print *, 'Enter output (PDF) filename :'
Read(*, 11) out

11 Format(a20)

! Open the input and output files
Open(1, File=in, Status='Old')
Open(2, File=out, Status='Unknown')

lbin = 100
i = 0

! Read data from input file until the end
do while (.true.)
    i = i + 1
    read(1, *, end=1000) a, xl(i)
end do 

1000 continue

niter = i - 1
xlmin = xl(1)
xlmax = xl(1)

! Find minimum and maximum values of xl
do i = 2, niter
    if (xl(i) .lt. xlmin) xlmin = xl(i)
    if (xl(i) .gt. xlmax) xlmax = xl(i)
end do

dn = (xlmax - xlmin) / dfloat(lbin)

! Process the data and write to output file
do i = 1, niter
    xx = (xl(i) - xlmin) / dn
    nxx = xx + 1
    ey(nxx) = ey(nxx) + 1
    xle = xle + xx
end do

sss = 0

! Find the maximum value in ey array
do i = 1, lbin
    if (ey(i) .gt. sss) sss = ey(i)
end do

do i = 1, lbin
    ! Write transformed data to output file
    write(2, *) xlmin + (i - 1) * dn, log(ey(i) / sss), ey(i) / sss, log(ey(i))
end do

end
