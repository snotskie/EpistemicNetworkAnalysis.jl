# TODO svd
# TODO ac-svd
# TODO mr1
# TODO mmr1
# TODO different plot options
# TODO saving plots
# TODO model fit testing
# TODO other (HLM, etc.) post-plot testing (with HypothesisTesting lib)

# using EpistemicNetworkAnalysis

function temp_example()
    RSdata = ena_dataset("RS.data")
    # display(first(RSdata, 6))
    # println()

    codes = [:Data,
             :Technical_Constraints,
             :Performance_Parameters,
             :Client_and_Consultant_Requests,
             :Design_Reasoning,
             :Collaboration]

    conversations = [:Condition, :GameHalf, :GroupName]
    units = [:Condition, :GameHalf, :UserName]
    # conversations = [:Condition, :GroupName]
    # units = [:Condition, :UserName]
    groupVar = :Condition
    controlGroup = "FirstGame"
    treatmentGroup = "SecondGame"
    # confounds = [:CONFIDENCE_Change]
    confounds = [:GameHalf]

    # myRotation = SVDRotation()
    # myRotation = MeansRotation(groupVar, controlGroup, treatmentGroup)
    myRotation = TwoGroupRotation(:Condition, "FirstGame", "SecondGame",
                                  :GameHalf, "First", "Second",
                                  [])

    myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=myRotation)
    display(myENA)

    # myArtist = MeansArtist(:Condition, "FirstGame", "SecondGame")
    # myArtist = MeansArtist(:GameHalf, "First", "Second")
    # myArtist = WindowsArtist(:Condition, "FirstGame", "SecondGame",
    #                          :GameHalf, "First", "Second")
    myArtist = TVRemoteArtist(:Condition, "FirstGame", "SecondGame",
                             :GameHalf, "First", "Second")
    scene = plot(myENA,
        showprojection=true,
        unitscale=0,
        artist=myArtist
    )
    display(scene)
end

# MR1 in Julia: X=MATCH, Y=APPROX
# │ Row │ relationship                                          │ thickness │ weight_x    │ weight_y    │
# │     │ Symbol                                                │ Float64   │ Float64     │ Float64     │
# ├─────┼───────────────────────────────────────────────────────┼───────────┼─────────────┼─────────────┤
# │ 1   │ Technical_Constraints_Performance_Parameters          │ 0.688247  │ -0.105832   │ 0.0671277   │
# │ 2   │ Data_Client_and_Consultant_Requests                   │ 0.179897  │ -0.00450161 │ 0.186166    │
# │ 3   │ Client_and_Consultant_Requests_Collaboration          │ 0.0589685 │ 0.00923557  │ 0.0601898   │
# │ 4   │ Design_Reasoning_Collaboration                        │ 0.223918  │ -0.150854   │ 0.274862    │
# │ 5   │ Data_Performance_Parameters                           │ 0.605124  │ 0.211841    │ 0.438982    │
# ⋮
# │ 10  │ Performance_Parameters_Client_and_Consultant_Requests │ 0.200761  │ 0.0701007   │ 0.219654    │
# │ 11  │ Performance_Parameters_Design_Reasoning               │ 0.740363  │ 0.532729    │ 0.103073    │
# │ 12  │ Data_Technical_Constraints                            │ 0.748241  │ -0.397539   │ -0.00707982 │
# │ 13  │ Technical_Constraints_Client_and_Consultant_Requests  │ 0.188269  │ -0.192893   │ 0.111187    │
# │ 14  │ Performance_Parameters_Collaboration                  │ 0.148319  │ -0.170993   │ 0.148614    │
# │ 15  │ Technical_Constraints_Design_Reasoning                │ 1.0       │ -0.120342   │ -0.698476   │

# MR1 in R: X=MATCH, Y=APPROX
# codes	MMR1	SVD2	SVD3	SVD4	SVD5	SVD6	SVD7	SVD8	SVD9	SVD10	SVD11	SVD12	SVD13	SVD14	SVD15
# <ena.mtdt>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>
# Data & Technical.Constraints                           	-0.397538600	0.05923297	0.40699051	-0.005182994	-0.338915642	-0.472034727	0.32207047	0.21063367	-0.29269678	0.099797253	-0.08582129	-0.16219275	-0.054639072	-0.21985444	-0.07991101
# Data & Performance.Parameters                          	0.211841163	-0.37998207	0.37804571	0.233328193	-0.044258797	0.383456650	-0.25242057	0.18217548	-0.20514222	0.544424187	-0.09800460	-0.08652517	-0.029313420	-0.09940494	-0.03407776
# Technical.Constraints & Performance.Parameters         	-0.105831538	-0.05668919	-0.14853470	0.627294404	0.003349060	0.049316434	-0.12031951	0.24644256	0.05883236	-0.375837464	0.10140539	0.01578707	0.286106882	-0.48625025	-0.13721283
# Data & Client.and.Consultant.Requests                  	-0.004501614	-0.22843294	-0.28790271	0.009198136	-0.283118142	0.010855556	0.43714985	-0.20046371	0.36658204	0.451023653	0.30249858	0.06131419	0.203979845	-0.13266755	-0.24695775
# Technical.Constraints & Client.and.Consultant.Requests	-0.192892945	-0.15848496	-0.22193186	-0.113598050	-0.226016062	-0.077663415	-0.28190795	-0.26316115	0.07227580	0.067262393	-0.52941286	-0.32831374	0.435371728	-0.07744279	0.27624983
# Performance.Parameters & Client.and.Consultant.Requests	0.070100706	-0.25882301	-0.26538405	0.076935825	-0.366915296	0.110423309	0.13585896	0.10344656	0.14722925	-0.191649266	-0.40787989	0.07875582	-0.660558139	-0.08902350	0.02308940
# Data & Design.Reasoning                                	0.210728648	0.03703219	0.39195592	-0.487512208	-0.137660953	0.006622412	-0.22421129	0.09612026	0.44402285	-0.167368170	-0.07053463	0.19245996	0.100424907	-0.39984655	-0.21742464
# Technical.Constraints & Design.Reasoning               	-0.120342195	0.64918067	-0.30603967	-0.127998170	0.024509942	0.290335457	-0.09437737	0.06544680	-0.11628437	0.335516547	-0.05245064	-0.12884213	-0.186278658	-0.41347327	-0.06962728
# Performance.Parameters & Design.Reasoning              	0.532729141	-0.13575026	-0.16617102	-0.018855325	0.427843316	-0.504853911	0.10212355	-0.05147660	-0.10352927	0.109447160	-0.12206028	-0.22065635	-0.094377832	-0.35037090	-0.02906697
# Client.and.Consultant.Requests & Design.Reasoning      	0.066135237	-0.21682875	-0.36787245	-0.294782734	-0.334587696	-0.198334678	-0.40289525	0.19276858	-0.44147600	-0.006602485	0.36044813	0.22151902	0.016747268	-0.03346495	-0.01661288
# Data & Collaboration                                   	-0.256853247	-0.21779487	0.13654506	-0.052971867	0.058049164	0.041269363	-0.21304006	-0.36787171	0.12930740	-0.112355497	0.47550460	-0.43690614	-0.374891947	-0.24271172	0.18320478
# Technical.Constraints & Collaboration                  	-0.527850970	-0.22050719	-0.13053958	-0.048746597	0.447470595	-0.224032089	-0.25632124	0.25374205	0.27256341	0.259480030	-0.12688810	0.21268903	-0.149641545	0.07274207	-0.20539058
# Performance.Parameters & Collaboration                 	-0.170992899	-0.13867524	0.06973731	-0.032378426	0.144294581	0.138583776	0.06174489	-0.60775920	-0.40569130	-0.102491841	-0.19098974	0.41277976	-0.008310563	-0.23129299	-0.31127022
# Client.and.Consultant.Requests & Collaboration         	0.009235566	-0.08319042	-0.09352695	-0.148545480	0.001183287	0.147285482	-0.01388922	0.08100281	-0.10111698	-0.202307248	-0.05644526	-0.53520580	0.036215709	0.24600858	-0.72877723
# Design.Reasoning & Collaboration                       	-0.150853866	-0.30682848	-0.08531594	-0.404381180	0.279871548	0.371516733	0.42664313	0.33554681	-0.15931847	-0.154015690	0.01346263	-0.06203900	0.168863909	-0.20220285	0.28070155

# MMR1 in Julia: X=MATCH, Y=APPROX
# │ Row │ relationship                                          │ thickness  │ weight_x   │ weight_y │
# │     │ Symbol                                                │ Real       │ Float64    │ Real     │
# ├─────┼───────────────────────────────────────────────────────┼────────────┼────────────┼──────────┤
# │ 1   │ Technical_Constraints_Performance_Parameters          │ -0.240922  │ -0.465829  │ 0        │
# │ 2   │ Data_Client_and_Consultant_Requests                   │ -0.0949388 │ -0.183566  │ 0        │
# │ 3   │ Client_and_Consultant_Requests_Collaboration          │ 0.0184638  │ 0.0357002  │ 0        │
# │ 4   │ Design_Reasoning_Collaboration                        │ -0.112032  │ -0.216616  │ 0        │
# │ 5   │ Data_Performance_Parameters                           │ 0.027958   │ 0.0540576  │ 0        │
# │ 6   │ Data_Collaboration                                    │ -0.0306977 │ -0.0593549 │ 0        │
# │ 7   │ Data_Design_Reasoning                                 │ 0.210225   │ 0.406475   │ 0        │
# │ 8   │ Client_and_Consultant_Requests_Design_Reasoning       │ 0.00495638 │ 0.00958328 │ 0        │
# │ 9   │ Technical_Constraints_Collaboration                   │ -0.18433   │ -0.356407  │ 0        │
# │ 10  │ Performance_Parameters_Client_and_Consultant_Requests │ -0.0984026 │ -0.190264  │ 0        │
# │ 11  │ Performance_Parameters_Design_Reasoning               │ -0.128562  │ -0.248578  │ 0        │
# │ 12  │ Data_Technical_Constraints                            │ 0.206711   │ 0.399682   │ 0        │
# │ 13  │ Technical_Constraints_Client_and_Consultant_Requests  │ 0.00118493 │ 0.00229109 │ 0        │
# │ 14  │ Performance_Parameters_Collaboration                  │ -0.0143781 │ -0.0278005 │ 0        │
# │ 15  │ Technical_Constraints_Design_Reasoning                │ 0.19613    │ 0.379223   │ 0        │

# MMR1 in R: X=MATCH, Y=APPROX
# codes	MMR1	SVD2	SVD3	SVD4	SVD5	SVD6	SVD7	SVD8	SVD9	SVD10	SVD11	SVD12	SVD13	SVD14	SVD15
# <ena.mtdt>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>	<en.dmnsn>
# Data & Technical.Constraints                           	-0.399681667	0.296224528	0.327010445	-0.1567488049	-0.21172555	-0.566148677	-0.083361574	0.21705996	-0.35285654	-0.09118961	-0.102869944	-0.06294312	-0.02668616	-0.22967246	0.027466941
# Data & Performance.Parameters                          	-0.054057571	-0.305790017	0.465749254	-0.1993600801	-0.22479115	0.569337486	-0.061908898	0.31989722	-0.26484499	0.29616906	-0.017034523	0.02144203	0.01119761	-0.08444631	-0.013956302
# Technical.Constraints & Performance.Parameters         	0.465828927	0.129201911	-0.144311989	-0.3654271273	-0.31497380	0.098176603	0.186238015	-0.20986049	-0.14467891	-0.22952722	-0.007185742	-0.07518060	0.27216429	-0.49156383	0.171310489
# Data & Client.and.Consultant.Requests                  	0.183566469	0.063765859	0.055234846	0.2902397885	-0.26927129	-0.257923674	0.132613175	0.24391549	0.34276208	0.56535951	0.287806466	-0.12793479	0.19071908	-0.12559359	0.265273388
# Technical.Constraints & Client.and.Consultant.Requests	-0.002291093	0.176718866	0.092665123	0.3414104748	-0.20868506	0.087271463	-0.327911998	-0.19646909	0.09692069	0.09158468	-0.519160332	0.08811946	0.50503708	-0.03103826	-0.308586175
# Performance.Parameters & Client.and.Consultant.Requests	0.190263854	-0.027195407	0.065268957	0.2613759976	-0.38013940	-0.033548642	0.304296212	-0.01777236	0.08574252	-0.03306255	-0.389373204	0.27356123	-0.63614110	-0.09742152	-0.053670366
# Data & Design.Reasoning                                	-0.406474782	-0.282835769	0.225532679	0.1000476404	0.23716831	-0.010404305	0.273741140	-0.53033389	0.08543952	0.13014543	-0.010120805	0.17717711	0.09013327	-0.40873020	0.224602031
# Technical.Constraints & Design.Reasoning               	-0.379222994	0.304494818	-0.577259337	0.1352366357	0.06136084	0.341306209	-0.004671212	0.23485497	-0.09056925	0.18860986	-0.069996640	-0.05863408	-0.13645996	-0.40704379	0.010092570
# Performance.Parameters & Design.Reasoning              	0.248577848	-0.598031926	-0.200712883	0.0761017641	0.22514921	-0.309853524	-0.370710451	0.19210327	-0.13038549	0.05639034	-0.237866514	-0.10763154	-0.05781486	-0.34894393	-0.002958485
# Client.and.Consultant.Requests & Design.Reasoning      	-0.009583280	-0.109508110	0.009078111	0.6043041089	-0.26197980	0.090060807	-0.177429102	-0.13297828	-0.38408697	-0.31765849	0.493616674	-0.02077372	-0.01628514	-0.06160293	0.012758521
# Data & Collaboration                                   	0.059354859	0.167906304	0.301114793	-0.0055948873	0.08373052	0.122983577	-0.235060557	-0.22377748	0.30114211	-0.02212033	0.072020606	-0.64994037	-0.34852845	-0.25372125	-0.213273487
# Technical.Constraints & Collaboration                  	0.356407044	0.393911493	0.131101576	0.0447008046	0.34561275	0.004294033	-0.238654784	-0.24053809	-0.40948620	0.40070168	0.009443915	0.23439163	-0.19456805	0.02295204	0.206610142
# Performance.Parameters & Collaboration                 	0.027800478	0.146833655	0.197310615	0.0002641221	0.12017065	0.122446785	-0.380355164	0.27927314	0.43879411	-0.36677014	0.094320242	0.44421410	-0.02296753	-0.24198410	0.307035564
# Client.and.Consultant.Requests & Collaboration         	-0.035700242	0.009964151	0.053865168	0.2037269953	0.03886470	0.134481309	0.046252764	0.05355123	-0.07852313	-0.18350720	-0.404360896	-0.40731266	0.04267076	0.24485034	0.706110992
# Design.Reasoning & Collaboration                       	0.216616149	0.136563369	0.255472148	0.3003715872	0.47235614	0.038538893	0.483676052	0.35819640	-0.11440549	-0.19060082	-0.046737601	-0.05933565	0.19221642	-0.17477194	-0.263573960