# TODO svd
# TODO ac-svd
# TODO mr1
# TODO mmr1
# TODO different plot options

function temp_example()
    RSdata = ena_dataset("RS.data")
    print(first(RSdata, 6))

    codes = [:Technical_Constraints, :Performance_Parameters, :Client_and_Consultant_Requests, :Design_Reasoning, :Collaboration]
    conversations = [:GroupName, :ActivityNumber]
    unitVar = :UserName

    myENA = ENAModel(RSdata, codes, conversations, unitVar)
    print(myENA.unitModel)
    print(myENA.networkModel)
end