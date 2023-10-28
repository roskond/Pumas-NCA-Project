using PharmaDatasets

df = dataset("iv_sd_1")

using NCA

#Mapping df with read_nca function
ncapop = read_nca(df; 
                id = :id,
                time = :time,
                observations = :conc,
                amt = :amt,
                route = :route,
)

#There are 30 subjects in the dataset

using NCAUtilities
using CairoMakie

#=Generating an Observation vs Time Plot 
and log-scale equivalent of Subject 27=#

#fig = Figure()
#observations_vs_time(
#    fig[1,1],
#    ncapop[27];
#    axis = (;
#        xlabel = "Concentration (mg/L)",
#        ylabel = "Time (Hours)",
#    ),
#)
#observations_vs_time(
#    fig[1,2],
#    ncapop[27];
#    axis = (;
#        xlabel = "Concentration (mg/L))",
#        ylabel = "Time (Hours)",
#        yscale = log10,
#    ),
#)

fig1 = Figure()
observations_vs_time(
    ncapop[27];
    axis = (;
        xlabel = "Concentration (mg/L)",
        ylabel = "Time (Hours)",
    ),
)

fig2 = Figure()
observations_vs_time(
    ncapop[27];
    axis = (;
        xlabel = "Concentration (mg/L)",
        ylabel = "Time (Hours)",
        yscale = log10,
    ),
)

#Summary Plot
fig3 = Figure()
summary_observations_vs_time(
    ncapop;
    axis = (;
        xlabel = "Concentation (mg/L)",
        ylabel = "Time (Hours)",
        ),
)

#Subject-fit plots

fig4 = subject_fits(
    ncapop;
    axis = (; yscale = log10),
    rows = 5,
    columns = 6,
)

using NCA.Unitful
using Dates

#NCA Analysis
report_1 = run_nca(
    ncapop;
    sigdigits = 3,
    studyid = "STUDY-002",
    studytitle = "Pharmacokinetics of Drug Z",
    author = [("Roshan Kondapalli", "Formula-Y")],
    sponsor = "Pumas-AI",
    date = Dates.now(),
    conclabel = "Concentration (mg/L)",
#    grouplabel = "Dosage (mg/L)",
    timelabel = "Time (Hours)",
    versionnumber = v"0.1",
)

#Parameter Summary
parm = [:cmax, :auclast,:kel, :half_life]
parm_sum = summarize(report_1; parameters = parm)

#Parameter Distribution
Cmax = parameters_dist(report_1; parameter = :cmax)
AUClast = parameters_dist(report_1; parameter = :auclast)
Kel = parameters_dist(report_1; parameter = :kel)
Half_life = parameters_dist(report_1; parameter = :half_life)

#Subject Fit Plots
individual_fit = subject_fits(
    report_1;
    axis = (;
        xlabel = "Time (Hours)",
        ylabel = "Concentration (mg/L)",
        yscale = log10
        ),
    seperate = true,
    paginate = true,
    limit = 15,
    columns = 5,
    rows = 3,
    facet = (; combinelabels = true),
)

#Computing Individual parameters

#Terminal Elimination Rate Constant
λz = NCA.lambdaz(ncapop; threshold = 3)

#Coefficient of Determination 
λzr² = NCA.lambdazr2(ncapop)    

#Correlation Coefficient
λzr = NCA.lambdazr(ncapop)

# y-intercept of λz
λzy = NCA.lambdazintercept(ncapop)

#=AUC calculation in 
time intervals of 0 to 12 
and 12 to 24=#
AUCs = NCA.auc(ncapop; interval =
     [(0,12), (12,24)]) 

final_report = innerjoin(
    report_1.reportdf, AUCs, λz, λzr, λzr²,
    λzy; on = [:id], makeunique = true)

report(report_1, output = "nca_report")