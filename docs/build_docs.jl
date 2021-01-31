using Pkg
Pkg.update("EpistemicNetworkAnalysis")

using Weave
weave("GettingStarted.jmd"; doctype="md2html", out_path=:pwd)
weave("ENAModel.jmd"; doctype="md2html", out_path=:pwd)
weave("Rotations.jmd"; doctype="md2html", out_path=:pwd)
weave("Plotting.jmd"; doctype="md2html", out_path=:pwd)