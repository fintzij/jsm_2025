# calculate restricted mean progression free survival time
function calc_rmpfst(paths)
    # get event times
    died = map(x -> x.states[2] != 1, paths)
    gaptimes = diff([0.0; sort(map(x -> x.times[2], paths[findall(died)])); maximum(map(x -> x.times[2], paths))])

    # counters
    rmpfst = 0.0
    nsurv  = length(paths)

    for k in 1:(length(gaptimes)-1)
        rmpfst += (nsurv - 0.5) * gaptimes[k] 
        nsurv -= 1
    end

    # add final increment
    rmpfst += nsurv * last(gaptimes)

    # return restricted mean PFS time
    return rmpfst / length(paths)
end

# summarize paths
function summarize_paths(paths; roundests = false)

    # progression-free survival
    # percent who progress
    # percent who die after progression
    # percent who die without progressing
    pfs = mean(map(x -> all(x.states .== 1), paths))
    prog = mean(map(x -> 2 ∈ x.states, paths))
    die_wprog = mean(map(x -> all([2,3] .∈ Ref(x.states)), paths))
    die_noprog = mean(map(x -> (3 ∈ x.states) & !(2 ∈ x.states), paths))

    # restricted mean progression-free survival time
    rmpfst = calc_rmpfst(paths)

    ests = roundests ? 
        (prog_free_surv = round(pfs, sigdigits = 2),
            progress = round(prog, sigdigits = 2),
            die_wprog = round(die_wprog, sigdigits = 2),
            die_noprog = round(die_noprog, sigdigits = 2),
            rmpfst = round(rmpfst, sigdigits = 2)) :
        (prog_free_surv = pfs, 
            progress = prog,
            die_wprog = die_wprog, 
            die_noprog = die_noprog, 
            rmpfst = rmpfst)

    return ests
end

# helper to create observation times for each subj
function make_obstimes()    
    # observation times
    collect(0:0.25:1) .+ vcat([0.0; (rand(Beta(5,5), 4) * 0.2 .- 0.1)])
end

# make data
# observation times for each subject
function make_data(nsubj, ntimes)
    # observation times for each subject
    obstimes = [make_obstimes() for i in 1:nsubj]

    # dataset - sets up simulation times and initial state
    data = DataFrame(
    id = repeat(collect(1:nsubj), inner = ntimes),
    tstart = reduce(vcat, map(x -> x[Not(end)], obstimes)),
    tstop = reduce(vcat, map(x -> x[Not(1)], obstimes)),
    statefrom = fill(1, nsubj * ntimes),
    stateto = fill(1, nsubj * ntimes),
    obstype = fill(2, nsubj * ntimes))

  return data
end

# get incidence and prevalence
function get_trajectory(paths)
    states = reduce(vcat, collect(map(x -> x.states[Not(1)], paths)))
    times = reduce(vcat, collect(map(x -> x.times[Not(1)], paths)))
    events = reduce(vcat, collect(map(x -> diff(x.states), paths)))
    ord = sortperm(times)

    states = states[ord]
    times = times[ord]
    events = events[ord]

    prev = DataFrame(time = [0;times],
                    healthy = nsubj,
                    progressed = 0,
                    dead = 0)

    inc = DataFrame(time = [0;times],
                    progression = 0,
                    healthy2death = 0,
                    progressed2death = 0)

    for k in 2:nrow(prev)
        if (states[k-1] == 2) & (events[k-1] != 0)
            prev.healthy[k]      = prev.healthy[k-1] - 1
            prev.progressed[k]   = prev.progressed[k-1] + 1
            prev.dead[k]         = prev.dead[k-1]
            inc.progression[k]   = inc.progression[k-1] + 1
            inc.healthy2death[k] = inc.healthy2death[k-1]
            inc.progressed2death[k] = inc.progressed2death[k-1]

        elseif events[k-1] == 2
            prev.healthy[k]      = prev.healthy[k-1] -1
            prev.progressed[k]   = prev.progressed[k-1]
            prev.dead[k]         = prev.dead[k-1] + 1
            inc.progression[k]   = inc.progression[k-1]
            inc.healthy2death[k] = inc.healthy2death[k-1] + 1
            inc.progressed2death[k] = inc.progressed2death[k-1]

        elseif events[k-1] == 1
            prev.healthy[k]      = prev.healthy[k-1]
            prev.progressed[k]   = prev.progressed[k-1] - 1
            prev.dead[k]         = prev.dead[k-1] + 1
            inc.progression[k]   = inc.progression[k-1]
            inc.healthy2death[k] = inc.healthy2death[k-1]
            inc.progressed2death[k] = inc.progressed2death[k-1] + 1
        elseif events[k-1] == 0
            prev.healthy[k]      = prev.healthy[k-1] 
            prev.progressed[k]   = prev.progressed[k-1]
            prev.dead[k]         = prev.dead[k-1] 
            inc.progression[k]   = inc.progression[k-1]
            inc.healthy2death[k] = inc.healthy2death[k-1] 
            inc.progressed2death[k] = inc.progressed2death[k-1]
        end
    end

    return prev, inc
end

function get_ests(model)
    paths = simulate(model; nsim = 1000, paths = true, data = false)
    summarize_paths(paths)
end

# asymptotic bootstrap
function asymptotic_bootstrap(model, pars, vcov, sims_per_subj, nboot)

    # draw parameters
    npars = length(pars)
    pardraws = zeros(Float64, npars)

    # SVD
    U = zeros(npars, npars) 
    D = zeros(npars)
    U,D = svd(vcov)

    # replace small negative singular values with zeros
    D[findall(D .< 0)] .= 0.0

    # matrix square root
    S = U * diagm(sqrt.(D))
    
    # initialize matrix of estimates
    ests = zeros(Float64, 5, nboot)

    # simulate paths under each set of parameters
    for k in 1:nboot
        # draw parameters
        pardraws[1:npars] = flatview(pars) .+ S * randn(npars)

        # set the parameters
        set_parameters!(model, VectorOfVectors(pardraws, model.parameters.elem_ptr))

        # simulate paths
        paths_sim = simulate(model; nsim = sims_per_subj, paths = true, data = false)

        # summarize paths
        ests[:,k] = collect(summarize_paths(paths_sim))
    end

    return mapslices(x -> quantile(x[findall(map(y -> (!ismissing(y) && (y != -1.0)), x))], [0.025, 0.975]), ests, dims = [2,])
end

function get_estimates(model)
    ests = get_ests(model)
    confints = asymptotic_bootstrap(model, flatview(model.parameters), model.vcov, 1, 1000) 
    return ests, confints
end