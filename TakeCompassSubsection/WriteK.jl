using JLD2, Polynomials
JLD2.@load "256K.jld2"

f = open("PERMX.txt","w") # do this once

write(f, "PERMX")
write(f, "\n")
ct = 0
for k = 1:256
    println("depth level = $(k)")
    for j = 1:256
        for i = 1:256
            global ct = ct + 1
            write(f, string(K[i,j,k]))
            write(f, " ")
            (ct%5 == 0) && write(f, "\n")
        end
    end
end

write(f, " /")
close(f)

f1 = open("PERMY.txt","w") # do this once

write(f1, "PERMY")
write(f1, "\n")
ct = 0
for k = 1:256
    println("depth level = $(k)")
    for j = 1:256
        for i = 1:256
            global ct = ct + 1
            write(f1, string(K[i,j,k]))
            write(f1, " ")
            (ct%5 == 0) && write(f1, "\n")
        end
    end
end

write(f1, " /")
close(f1)


f2 = open("PERMZ.txt","w") # do this once

write(f2, "PERMZ")
write(f2, "\n")
ct = 0
for k = 1:256
    println("depth level = $(k)")
    for j = 1:256
        for i = 1:256
            global ct = ct + 1
            write(f2, string(K[i,j,k]))
            write(f2, " ")
            (ct%5 == 0) && write(f2, "\n")
        end
    end
end

write(f2, " /")
close(f2)

f3 = open("PORO.txt","w") # do this once

write(f3, "PORO")
write(f3, "\n")
ct = 0
for k = 1:256
    println("depth level = $(k)")
    for j = 1:256
        for i = 1:256
            global ct = ct + 1
            p = Polynomial([-0.0314^2*K[i,j,k],2*0.0314^2*K[i,j,k],-0.0314^2*K[i,j,k],1.527^2])
            phi = minimum(real(roots(p)[findall(real(roots(p)).== roots(p))]))
            write(f3, string(phi))
            write(f3, " ")
            (ct%5 == 0) && write(f3, "\n")
        end
    end
end
write(f3, " /")
close(f3)