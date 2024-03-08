function Jac = Jacobian(f, x)

for f_i=1:size(f, 1)
    for x_i=1:size(x, 1)
        Jac(f_i, x_i) = diff(f(f_i), x(x_i));
    end
end

end

