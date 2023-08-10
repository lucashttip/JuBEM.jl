using Optim

function findmind_optim(points, source_point, csis_cont)

    lower = [-1.0,-1.0]
    upper = [1.0,1.0]
    res = optimize(Optim.only_fg!((F,G,x) -> fg!(F,G,x;points=points, source=source_point, csis_cont=csis_cont)), lower, upper,[0.0,0.0], Fminbox(BFGS()))

    return Optim.minimizer(res), minimum(res)
end