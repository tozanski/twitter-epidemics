module TwitterPlots

using Twitter
using PyPlot
using KernelDensity
using TomekUtils

plot_dir = "plots"
function set_plot_dir(new_plot_dir::String)
  global plot_dir
  plot_dir = new_plot_dir
end


function degree_ecdf_plot(adj, datasetname::String; fig=figure())
  deg_out = vec(sum(adj,dims=2))
  deg_in = vec(sum(adj,dims=1))

  x_out,p_out = arrecdf(deg_out)
  x_in, p_in = arrecdf(deg_in)

#  fig = figure()
  semilogx(x_in, p_in, label="in-degree")
  semilogx(x_out, p_out, label="out-degree")
  #title("Data set \"$datasetname\"")
  legend(loc="lower right")
  xlabel("degree")
  ylabel("cdf")
  xlim(1, exp10(ceil(log10(max(x_in[end], x_out[end])))))
  ylim(0,1)
  savefig(plot_dir*"/degree_ecdf_$(datasetname).pgf")
  close(fig)
  fig
end

function inoutdeg_corr_plot(adj,adj2,datasetname::String, secondname::String; fig=figure())
  loglog( vec(sum(adj,dims=1)), vec(sum(adj,dims=2)),".", markersize=1.0, label="original")
  loglog( vec(sum(adj2,dims=1)), vec(sum(adj2,dims=2)),".", markersize=1.0, label=secondname)
  #title("Data set \"$(datasetname)\", original vs $secondname")
  legend()
  xlabel("in-degree")
  ylabel("out-degree")
  savefig(plot_dir*"/inoutdeg_corr_$(datasetname)_$(secondname).png")
  close(fig)
  fig
end

function inoutdeg_corr_plot(adj, datasetname::String; fig=figure())
  loglog( vec(sum(adj,dims=1)), vec(sum(adj,dims=2)),".", markersize=1.0, label="original")
  #loglog( vec(sum(adj2,dims=1)), vec(sum(adj2,dims=2)),".", markersize=1.0, label=secondname)
  #title("Data set \"$(datasetname)\", original vs $secondname")
  legend()
  xlabel("in-degree")
  ylabel("out-degree")
#  savefig(plot_dir*"/inoutdeg_corr_$(datasetname).png")
  close(fig)
  fig
end

function inoutdeg_corr_contour(adj, datasetname::String; fig=figure())
  log_deg_in  = sum(adj, dims=1) |> vec .|> log10 .|> x-> max(-1,x)
  log_deg_out = sum(adj, dims=2) |> vec .|> log10 .|> x-> max(-1,x)

  density = kde(
      (log_deg_in, log_deg_out),
      npoints=(2048,2048),
      boundary=((-3,7), (-3,7)),
      bandwidth=(0.1,0.1)
      )

  levels = [0.001; 0.002; 0.005; 0.01; 0.02; 0.05:0.05:1.0]
  #levels = [0.05, 0.25, 0.5, 1.0]
  cx = contour(density.x, density.y, density.density, levels=levels)
  axis("equal")
  xmax = 4
  ymax = 4
  xlim(-1.5, xmax)
  ylim(-1.5, ymax)
  xticks( 0:xmax, 0:floor(Int64, xmax) .|> x-> raw"$10^"*string(x)*raw"$")
  yticks( 0:ymax, 0:floor(Int64, ymax) .|> x-> raw"$10^"*string(x)*raw"$")
  #xlabel("in-degee")
  #ylabel("out-degree")
  #colorbar()
  savefig(plot_dir*"/inoutdeg_corr_contour_$datasetname.pgf")
  #fig
end

function gamma_plot(gammas, likelihoods, datasetname; fig=figure())
  semilogx(gammas, likelihoods)
  xlabel(raw"$\gamma$")
  ylabel("log-likelihood")
  #title("Data set \"$(datasetname)\"")
  savefig(plot_dir*"/gamma_likelihood_$(datasetname).pgf")
  close(fig)
  fig
end

function forward_epidemics_plot(data; fig=figure())
  transprob = data["transprob"]
  datasetname = data["datasetname"]

  line_est_th, = loglog( transprob, data["outbreak_theo"], "-", color = "C1")
  line_shu_th, = loglog( transprob, data["outbreak_theo_ind"], "--", color = "C2")
  line_sca_th, = loglog( transprob, data["outbreak_theo_sca"], ":", color = "C3")

  line_est, = loglog( transprob, data["outbreak_est"], "v", color = "C1")
  line_shu, = loglog( transprob, data["outbreak_shu"], "s", color = "C2")
  line_sca, = loglog( transprob, data["outbreak_sca"], "o", color = "C3")

  line_orig, = loglog( transprob, data["outbreak"], "^", color = "C0")


  #title("forward epidemics '$datasetname'")
  xlabel("transmission probability")
  ylabel("prob. to reach large size")

  lines = [
          line_orig,
          line_est, line_est_th,
          line_shu, line_shu_th,
          line_sca, line_sca_th]

  legend(lines,[
        "original",
        "2D", "theory",
        "product", "theory",
        "scalar", "theory"],
        loc = "lower right"
  )

  xlim(10^-3, 1)
  ylim(10^-3, 1)

  savefig(plot_dir*"/fwd_$datasetname.pgf")
  close(fig)
  fig
end

function backward_epidemics_plot(data; fig=figure())
  transprob = data["transprob"]
  datasetname = data["datasetname"]

  line_est_th, = loglog( transprob, data["prevalence_theo"], "-", color = "C1")
  line_shu_th, = loglog( transprob, data["prevalence_theo_ind"], "--", color = "C2")
  line_sca_th, = loglog( transprob, data["prevalence_theo_sca"], ":", color = "C3")

  line_est, = loglog( transprob, data["prevalence_est"], "v", color = "C1")
  line_shu, = loglog( transprob, data["prevalence_shu"], "s", color = "C2")
  line_sca, = loglog( transprob, data["prevalence_sca"], "o", color = "C3")

  line_orig, = loglog( transprob, data["prevalence"], "^", color = "C0")

  #title("backward epidemics '$datasetname'")
  xlabel("transmission probability")
  ylabel("prob. to reach large size")
  #legend(loc="lower right")
  lines = [
          line_orig,
          line_est, line_est_th,
          line_shu, line_shu_th,
          line_sca, line_sca_th]

  legend(lines,[
        "original",
        "2D", "theory",
        "product", "theory",
        "scalar", "theory"],
        loc = "lower right"
  )

  xlim(10^-3, 1)
  ylim(10^-3, 1)
  savefig(plot_dir*"/bwd_$datasetname.pgf")
  close(fig)
  fig
end


function scc_plot(data; fig=figure())
  transprob = data["transprob"]
  datasetname = data["datasetname"]

  line_est_th, = loglog( transprob, data["scc_theo"], "-", color = "C1")
  line_shu_th, = loglog( transprob, data["scc_theo_ind"], "--", color = "C2")
  line_sca_th, = loglog( transprob, data["scc_theo_sca"], ":", color = "C3")

  line_est, = loglog( transprob, data["scc_est"], "v", color = "C1")
  line_shu, = loglog( transprob, data["scc_shu"], "s", color = "C2")
  line_sca, = loglog( transprob, data["scc_sca"], "o", color = "C3")

  line_orig, = loglog( transprob, data["scc"], "^", color = "C0")

  #title("strongly connected component '$datasetname'")
  xlabel("transmission probability")
  ylabel("nodes belonging to largest SCC")
  #legend(loc="lower right")
  lines = [
          line_orig,
          line_est, line_est_th,
          line_shu, line_shu_th,
          line_sca, line_sca_th]

  legend(lines,[
        "original",
        "2D", "theory",
        "product", "theory",
        "scalar", "theory"],
        loc = "lower right"
  )
  xlim(10^-3, 1)
  ylim(10^-3, 1)
  savefig(plot_dir*"/scc_$datasetname.pgf")
  close(fig)
  fig
end

function resistances_plot(data; fig=figure())
  transprob = data["transprob"]
  datasetname = data["datasetname"]



  semilogx( transprob, data["resistance"], "^", label="original")
  semilogx( transprob, data["resistance_est"], "v", label="resampled")
  semilogx( transprob, data["resistance_shu"], "s", label="resampled shuffled")
  semilogx( transprob, data["resistance_sca"], "o", label="resampled scalar")

  semilogx( transprob, data["resistance_theo"], "-", label="fixed point equation")
  semilogx( transprob, data["resistance_theo_ind"], "--", label="fixed point equation")
  semilogx( transprob, data["resistance_theo_sca"], ":", label="fixed point equation")
  #title("resistance '$datasetname'")
  xlabel("transmission probability")
  ylabel("resistance")
  #legend(loc="upper left")
  ylim(10.0^-3, 1)
  savefig(plot_dir*"/resistance_$datasetname.pgf")
  close(fig)
  fig
end

function contacts_plot(data; fig=figure())
  transprob = data["transprob"]
  datasetname = data["datasetname"]

  line_est_th, = loglog( transprob, data["contact_theo"], "-", color = "C1")
  line_shu_th, = loglog( transprob, data["contact_theo_ind"], "--", color = "C2")
  line_sca_th, = loglog( transprob, data["contact_theo_sca"], ":", color = "C3")

  line_est, = loglog( transprob, data["contact_est"], "v", color = "C1")
  line_shu, = loglog( transprob, data["contact_shu"], "s", color = "C2")
  line_sca, = loglog( transprob, data["contact_sca"], "o", color = "C3")

  line_orig, = loglog( transprob, data["contact"], "^", color = "C0")

  #title("contact probability '$datasetname'")
  xlabel("transmission probability")
  ylabel("contact probability")
  #legend(loc="lower right")
  lines = [
          line_orig,
          line_est, line_est_th,
          line_shu, line_shu_th,
          line_sca, line_sca_th]

  legend(lines,[
        "original",
        "2D", "theory",
        "product", "theory",
        "scalar", "theory"],
        loc = "lower right"
  )
  xlim(10^-3, 1)
  ylim(10^-3, 1)
  savefig(plot_dir*"/contact_$datasetname.pgf")
  close(fig)
  fig
end

function original_epidemics_plot(data; fig=figure())
  transprob = data["transprob"]

  line_scc_th, = loglog( transprob, data["scc_theo"], "-", color = "C0")
  line_out_th, = loglog( transprob, data["outbreak_theo"], "--", color = "C1")
  line_pre_th, = loglog( transprob, data["prevalence_theo"], ":", color = "C2")
  line_con_th, = loglog( transprob, data["contact_theo"], "-.", color = "C3")

  line_scc, = loglog( transprob, data["scc"], "s", color="C0")
  line_out, = loglog( transprob, data["outbreak"], "^", color="C1")
  line_pre, = loglog( transprob, data["prevalence"], "v", color="C2")
  line_con, = loglog( transprob, data["contact"], "o", color="C3")

  #title("'$(data["datasetname"])' - original graph")
  xlabel("transmission probability")
  ylabel("probability")

  lines = [
          line_scc, line_scc_th,
          line_out, line_out_th,
          line_pre, line_pre_th,
          line_con, line_con_th]

  legend(lines,[
        "scc", "theory",
        "outbreak", "theory",
        "prevalence", "theory",
        "reach", "theory"],
        loc = "lower right"
  )
  xlim(10^-3, 1)
  ylim(10^-3, 1)
  savefig(plot_dir*"/orig_$(data["datasetname"]).pgf")
  close(fig)
  fig
end

function estimated_epidemics_plot(data; fig=figure())
  transprob = data["transprob"]
  loglog( transprob, data["scc_est"], "^", label="scc")
  loglog( transprob, data["outbreak_est"], "^", label="outbreak")
  loglog( transprob, data["prevalence_est"], "^", label="prevalence")
  loglog( transprob, data["contact_est"],"o", label="contact" )

  loglog( transprob, data["scc_theo"], "-", label="scc (theory)")
  loglog( transprob, data["outbreak_theo"], "--", label="outbreak (theory)")
  loglog( transprob, data["prevalence_theo"], ":", label="prevalence (theory)")
  loglog( transprob, data["contact_theo"], "-.", label = "contact (theory)")

  #title("'$(data["datasetname"])' - estimated graph ")
  xlabel("transmission probability")
  ylabel("nodes reaching giant size")
  #legend(loc="lower right")
  legend()
  xlim(10^-3, 1)
  ylim(10^-3, 1)
  savefig(plot_dir*"/estimated_$(data["datasetname"]).pgf")
  close(fig)
  fig
end

function estimated_independent_epidemics_plot(data; fig=figure())
  transprob = data["transprob"]
  loglog( transprob, data["scc_shu"], "^", label="scc")
  loglog( transprob, data["outbreak_shu"], "^", label="outbreak")
  loglog( transprob, data["prevalence_shu"], "^", label="prevalence")
  loglog( transprob, data["contact_shu"], "o", label="contact")

  loglog( transprob, data["scc_theo_ind"], "-", label="scc (theory)")
  loglog( transprob, data["outbreak_theo_ind"], "--", label="outbreak (theory)")
  loglog( transprob, data["prevalence_theo_ind"], ":", label="prevalence (theory)")
  loglog( transprob, data["contact_theo_ind"], "-.", label = "contact (theory)")

  #loglog( transprob, data["outbreak_shu"].*data["prevalence_shu"], "-", label="prevalence (prod theory)")


  #title("'$(data["datasetname"])' - estimated graph with independent edges")
  xlabel("transmission probability")
  ylabel("nodes reaching giant size")
  #legend(loc="upper left")
  legend()
  ylim(10.0^-3, 1)
  savefig(plot_dir*"/estimated_independent_$(data["datasetname"]).pgf")
  close(fig)
  fig
end

function estimated_scalar_epidemics_plot(data; fig=figure())
  transprob = data["transprob"]
  loglog( transprob, data["scc_sca"], "^", label="scc")
  loglog( transprob, data["outbreak_sca"], "^", label="outbreak")
  loglog( transprob, data["prevalence_sca"], "^", label="prevalence")
  loglog( transprob, data["contact_sca"], "o", label="contact")

  loglog( transprob, data["scc_theo_sca"], "-", label="scc (theory)")
  loglog( transprob, data["outbreak_theo_sca"], "--", label="outbreak (theory)")
  loglog( transprob, data["prevalence_theo_sca"], ":", label="prevalence (theory)")
  loglog( transprob, data["contact_theo_sca"], "-.", label = "contact (theory)")

  #title("'$(data["datasetname"])' - estimated graph with scalar kernel")
  xlabel("transmission probability")
  ylabel("nodes reaching giant size")
  #legend(loc="upper left")
  legend()
  ylim(10.0^-3, 1)
  savefig(plot_dir*"/estimated_scalar_$(data["datasetname"]).pgf")
  close(fig)
  fig
end

function outbreak_prevalence_dichotomy_plot(data; fig=figure())
  datasetname = data["datasetname"]
  figure(figsize=(20,20))
  plot(data["transprob"], data["outbreak"] - data["prevalence"], "o", label="original")
  plot(data["transprob"], data["outbreak_theo"] - data["prevalence_theo"], label="theory")

  #title("'$datasetname' - outbreak - prevalence dichotomy")
  xlabel("transmission probability")
  ylabel("outbreak - prevalence")
  legend()
  savefig(plot_dir*"/dichotomy_$datasetname.pgf")
end

end