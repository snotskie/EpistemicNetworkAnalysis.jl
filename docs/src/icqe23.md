# ICQE23

## What to Bring

For this workshop you'll need:

- 💻 Your laptop, setup as described below so we can hit the ground running
- 🔌 Laptop charger
- 🎉 A can-do attitude and an unshakeable curiousity
- 🗂️ Any data you would like to work on, help thinking about, etc.

## Before the Conference: Setup and Survey

Before the conference begins:

1. [Install Julia](https://julialang.org/downloads/)
2. [Install VS Code](https://code.visualstudio.com/download)
3. [Install Julia VS Code extension](https://code.visualstudio.com/docs/languages/julia#_getting-started)
4. Install `EpistemicNetworkAnalysis.jl` by:
    1. Create a new file in VS Code
    2. Save it with the name `setup.jl`
    3. Copy the following code into that file and save the changes:

            # Install ENA package
            using Pkg
            Pkg.add(url="https://github.com/snotskie/EpistemicNetworkAnalysis.jl")

    4. Run the file by pressing Run / Run without Debugging. Note, the first time you run this code may time some time to complete
5. Once everything is installed, test that it works by:
    1. Create a new file in VS Code
    2. Save it with the name `main.jl`
    3. Copy the following code into that file and save the changes:

            # Load ENA package
            using EpistemicNetworkAnalysis

            # Load sample dataset, codes from my first year on hormone replacement therapy
            data = loadExample("transitions")

            # Derive some new codes based on old ones
            deriveAnyCode!(data, :BODY, :Changes, :Mood, :Oily, :Dysphoria, :Cry)
            deriveAnyCode!(data, :REFLECT, :Identity, :Longing, :Dream, :Childhood, :Family, :Name, :Letter, :Doubt, :Religion)
            deriveAnyCode!(data, :LEARN, :WWW, :Experiment, :Recipe)
            deriveAnyCode!(data, :PROGRESS, :Strangers, :Passed, :Out, :Affirmation)

            # Add new columns for splitting the year's data in half, third, ...
            data[!, :All] .= "All"
            data[!, :Half] .= "First"
            data[183:end, :Half] .= "Second"
            data[!, :Third] .= "First"
            data[122:243, :Third] .= "Second"
            data[244:end, :Third] .= "Third"
            data[!, :Fourth] .= "First"
            data[92:183, :Fourth] .= "Second"
            data[184:275, :Fourth] .= "Third"
            data[276:end, :Fourth] .= "Fourth"

            # List columns to use as codes, convos, and units
            codes = [:DoseTracking, :SkippedDose, :Happy, :NonHappy, :Sweets, :BODY, :REFLECT, :LEARN, :PROGRESS]
            conversations = [:Date]
            units = [:Date]

            # Run the model and plot it
            model = ENAModel(data, codes, conversations, units)
            p = plot(model)
            display(p)

    4. Run the file by pressing Run / Run without debugging. Note, the first time you run this code may time some time to complete
    5. If an ENA plot appears, congrats! Everything is setup and ready to go for the workshop!
5. [Complete the prior knowledge survey](https://forms.gle/8QVEy6NxnQCY9tBz9) to help us plan and tell us a little about yourself

If you have any problems and questions during setup at all, [please let me know](mailto:mariah.knowles@wisc.edu?subject=ICQE23+Setup+Question) and we'll get you squared away!

## Intro

This is a 2-hour tutorial workshop

Workshop Goals:

- Give QE researchers the power to choose and develop custom rotations to fit their research aims, by introducing them to open-source ENA-based tools
- Learners will be able to (a) create and interpret ENA models that have rotations beyond SVD and means rotation
- (b) understand the steps of the ENA algorithm and the connections between rotation choice and research aims
- (c) use Github issues to get assistance, troubleshoot problems, and contribute to ENA API development

Workshop Audience:

- QE researchers
- already familiar with ENA
- a beginner's understanding of at least one programming language
- who want to use rotations beyond SVD and means rotation in their own research

Workshop structure:

- Live-coding / code-along model
- Three phases:
    1. Demo how to make several different models with different rotations with one running example
    2. Discussion on choosing the right rotation for the job, based on your gut sense of your research story
    3. Looking ahead discussion, sharing pathways for contributing and getting help after the workshop

Code of Conduct:

- Use welcoming and inclusive language
- Be respectful of different viewpoints and experiences
- Gracefully accept constructive criticism
- Focus on what is best for the community
- Show courtesy and respect towards other community members

[Full CoC](https://docs.carpentries.org/topic_folders/policies/code-of-conduct.html)

If you would like to report any misconduct, contact a Data Science Hub Facilitator: <facilitator@datascience.wisc.edu>

**Collaborative Notes**

We use Etherpad to write notes collaboratively, link to be given during the event

Notes and code will be added as the instructor presents the material

**Have conceptual questions?**

Get the instructor’s attention by raising your hand

**Encounter a bug?**

Get a helper’s attention by raising your hand or putting your name card on its side with the red facing up

Setup check: green side up if they completed setup before the workshop

Activity, introductions, keep running notes on the board:

- What are you excited to learn from this workshop?
- Or one question you have?

## Demo and Worked Example

**Dataset**

The dataset we'll be using are codes and metadata from my first year on hormone replacement therapy. I'm a transgender woman, Valentines Day 2020 I started my medical transition, and I kept a daily record in various ways on my phone that first year. Also that year I learned QE and started developing tools for ENA. As a test case, I coded my own data, [modeled it](https://link.springer.com/chapter/10.1007/978-3-030-93859-8_8), and more recently made it available. (Though, just the codes and metadata, not the text of the daily entries, given their private nature.)

**Getting Started**

We'll start where the setup instructions left off:

```julia
# Load ENA package
using EpistemicNetworkAnalysis

# Load sample dataset, codes from my first year on hormone replacement therapy
data = loadExample("transitions")

# Derive some new codes based on old ones
deriveAnyCode!(data, :BODY, :Changes, :Mood, :Oily, :Dysphoria, :Cry)
deriveAnyCode!(data, :REFLECT, :Identity, :Longing, :Dream, :Childhood, :Family, :Name, :Letter, :Doubt, :Religion)
deriveAnyCode!(data, :LEARN, :WWW, :Experiment, :Recipe)
deriveAnyCode!(data, :PROGRESS, :Strangers, :Passed, :Out, :Affirmation)

# Add new columns for splitting the year's data in half, third, ...
data[!, :All] .= "All"
data[!, :Half] .= "First"
data[183:end, :Half] .= "Second"
data[!, :Third] .= "First"
data[122:243, :Third] .= "Second"
data[244:end, :Third] .= "Third"
data[!, :Fourth] .= "First"
data[92:183, :Fourth] .= "Second"
data[184:275, :Fourth] .= "Third"
data[276:end, :Fourth] .= "Fourth"

# List columns to use as codes, convos, and units
codes = [:DoseTracking, :SkippedDose, :Happy, :NonHappy, :Sweets, :BODY, :REFLECT, :LEARN, :PROGRESS]
conversations = [:Date]
units = [:Date]

# Run the model and plot it
model = ENAModel(data, codes, conversations, units)
p = plot(model)
display(p)
```

**Conversation, Window Size, and Empty Units**

As is, each day is its own "conversation," meaning no connections are made from day to day. Instead, lets change the conversation to the whole year:

```julia
conversations = [:All] # :All is a column we made that just contains the word "All"
```

That's better, I guess. But let's inspect the model to see what's going on:

```julia
display(model)
```

Notice this part of the output:

```txt
 ModelConfig =
     (codes = [:DoseTracking,
               :SkippedDose,
               :Happy,
               :NonHappy,
               :Sweets,
               :BODY,
               :REFLECT,
               :LEARN,
               :PROGRESS],
      conversations = [:All],
      units = [:Date],
      unitFilter = ...,
      edgeFilter = ...,
      windowSize = Inf, <--- This part right here
      sphereNormalize = true,
      dropEmpty = false,
      recenterEmpty = false),
```

By default, the window size in this package is infinite, meaning any connection between any codes in the whole conversation are counted. Let's pick a more sensible window:

```julia
model = ENAModel(
    data, codes, conversations, units,
    windowSize=4
)
```

Much better. One more tweak. Let's tell the model what to do with empty units. David's convinced me that the best place to put them is in the mean or center of the plot. This has advantages for proper stats tests later

```julia
model = ENAModel(
    data, codes, conversations, units,
    windowSize=4,
    recenterEmpty=true
)
```

**Plot Tweaks**

This is data that moves over time. And other projects I know have data that moves over different continuous variables, like grades or so on. Let's color code our plot to make that pop:

```julia
p = plot(
    model,
    spectralColorBy=true
)
```

Alternatively, we could split the year's data in half and color code by that instead:

```julia
p = plot(
    model,
    groupBy=:Half # :Half is a column we added that says "First" then "Second"
)
```

Or we could split it in thirds (`groupBy=:Third`) or in fourths (`groupBy=:Fourth`)

Activity:

- Which version of the plot shows the most information?
- What "story" are you starting to see emerge about the data from that year?

**Rotations**

The first and second halves seem different. Let's rotate the model to make that the focus:

```julia
rotation = MeansRotation(:Half, "First", "Second")
model = ENAModel(
    data, codes, conversations, units,
    windowSize=4,
    recenterEmpty=true,
    rotateBy=rotation
)
```

Surely though, a person's experience over a year can't be summed up as just a before and an after, right? Let's split the year in thirds and make that three-way difference the focus:

```julia
rotation = MulticlassRotation(:Third)
```

Note the test statistic for the x- and y-axes:

```julia
display(model.embedding[1, :KruskalWallis_H])
display(model.embedding[2, :KruskalWallis_H])
```

If we *really* want to pull the groups apart and maximize that statistic, we can use a slightly different rotation:

```julia
rotation = LDARotation(:Third)
```

Note the drop in the variance explained. The plot is also hard to read, so let's zoom in:

```julia
p = plot(
    model,
    zoom=3
)
```

Activity:

- Which of those two models (Multiclass or LDA) do you prefer and why?

Let's try splitting the year in fourth next:

```julia
rotation = MulticlassRotation(:Fourth)
# or
rotation = LDARotation(:Fourth)
```

Wow, that's a lot of subplots. Let's pull out just one for now:

```julia
sp = plot(p.subplots[1], size=(600,600))
display(sp)
```

Activity:

- We've seen this data modeled a lot of different ways now. Which codes seem to really be driving the story the most?

Let's model that assumption directly:

```julia
rotation = TopicRotation("HRT?", [:SkippedDose], [:Happy])
# ...
p = plot(
    model,
    spectralColorBy=true
)
# or
p = plot(
    model,
    groupBy=:Fourth
)
```

That seems to be pretty close. The spectral plot still has a lot of red in the middle though. Being familiar with this data, and having lived it, I know that there's a little more to the story than just those two codes

```julia
rotation = TopicRotation("HRT", [:SkippedDose, :DoseTracking], [:Happy, :PROGRESS])
```

Much better

That was all a lot of work. What if we just wanted to ask, show me the model that focuses on time as a linear scale?

```julia
using Pkg
Pkg.add("GLM")
using GLM
# ...
rotation = FormulaRotation(
    LinearModel, @formula(y ~ 1 + Day), 2, nothing
)
```

## Github

`EpistemicNetworkAnalysis.jl` lives on [GitHub](https://github.com/snotskie/EpistemicNetworkAnalysis.jl)

But not just that package. All of Julia lives on GitHub. Every Julia package, and even the language itself, by design are hosted on GitHub. This can make things difficult for people living in [some countries where GitHub is occasionally blocked or slowed](https://en.wikipedia.org/wiki/Censorship_of_GitHub). By and far though, this allows you to:

- See and report issues for all Julia packages, as well as propose your own fixes
- See the entire history of development and discussion for all Julia packages
- Easily find documentation for most Julia packages
- Engage in discussion with other users of some of your favorite packages
- Directly ask developers for help

Julia has from the very beginning made open source development a part of how the language itself works. I really appreciate this because it builds community right into everything I do. You never code alone, and you never learn to code alone

While Julia takes this to the limit, other QE resources are open to some extent:

- rENA lets you see issues and development history on [gitlab](https://gitlab.com/epistemic-analytics/qe-packages/rENA/-/issues)
- Several workshops shared on [qesoc](https://qesoc.org/) via Google Drive:
    - ROCK
    - Intro to ENA
    - Advanced ENA
    - rENA
    - Intro to Automated Coding

Those workshops are *not* shared open source in the same sense as GitHub and gitlab, in that we can't provide feedback, suggest changes, ask questions, or review histories

Activity:

1. Summarize with your neighbor in 10 words or less: What do we care about when we talk about *open source* research software?
2. What other QE open resources am I missing? Other not-specifically-QE-but-still-useful open resources?

## Julia Basics

Variables in Julia are most similar to Python

- `i = 1.0` creates a new variable
- `print(i)` shows its value
    - so does `println(i)`
    - and `display(i)`
    - and `dump(i)`
    - and `show(i)`
    - and `@show i`, which also adds the name of the variable to the output
    - use `print`/`println` for text-only stuff and `display` for "interactive" stuff
- `typeof(i)` gives the data type of `i`
- `supertype(Float64)` gives the type above that type in the type hierarchy
- `supertypes(Float64)` gives all the types above it, in order
- `subtypes(Real)` gives all the immediate types under `Real`

The `TypeTree` package can visualize it all the way down

```julia
using Pkg
Pkg.add("TypeTree")
using TypeTree
print.(tt(Real))
```

We'll see this again later

A review of numbers and strings:

- `1` is an int
- `1.0` is a float
- `1 + 1im` is a complex
- `"hello"` is a string
    - so is `"""hello"""`
    - but not `'hello'`
- `'h'` is a single character, like in C and Java
- `:Hello` is a symbol

We'll use symbols for column names and strings for column values. If you have `:Hello` in five places in your code, they always resolve to the exact same thing in memory. It is the name for the exact same concrete label. But if you have `"hello"` in five places, those could resolve to five different things in memory, like having five people named `"David"` at the same meeting

Arrays in Julia work similar to Python, with one main difference:

- `a = [1, 2, 3]` creates an array
- `a[1]` is the *first* item in the array, like in R

Julia has four consoles in a single REPL:

- `;` brings up the shell console
- `?` brings up the help console
- `]` brings up the package console
- and backspacing takes you back to the main console

We'll use the main console and the help console the most

Julia has robust type and package systems. I've taught a [4-hour workshop](https://carpentries-incubator.github.io/julia-novice/) on that alone. A lot of what we're doing today benefits from those type and package systems working in the background. We'll just stick to the basics

- `using Pkg` loads the `Pkg` library into memory and gives us access to all the tools it exports
- `Pkg.add(...)` installs a package into Julia's central package repo on your computer. If it requires other packages, those will be downloaded too, but won't automatically be marked "added," meaning they exist only to fulfill that dependency and Julia is free to clean it up if it stops being required

`EpistemicNetworkAnalysis.jl` depends on `DataFrames` and `CSV`.

This will fail until we `Pkg.add(...)` them (assuming you haven't done that yourself already):

```julia
using DataFrames
using CSV
```

Activity, fill in the blanks:

```julia
[___] Pkg
Pkg.[___]("Statistics")
Pkg.[___]("HypothesisTests")
Pkg.[___]("Plots")

[___] Statistics
[___] HypothesisTests
[___] Plots
```

## EpistemicNetworkAnalysis.jl

ENA models in `EpistemicNetworkAnalysis.jl` have the following required and optional parameters:

```julia
ENAModel(
    # Required
    data::DataFrame,
    codes::Array{Symbol,1},
    conversations::Array{Symbol,1},
    units::Array{Symbol,1};

    # Optional
    rotation::AbstractLinearENARotation=SVDRotation(),
    unitFilter::Function=unit->true,
    edgeFilter::Function=edge->edge.kind == :undirected,
    windowSize::Real=Inf,
    sphereNormalize::Bool=true,
    dropEmpty::Bool=false,
    recenterEmpty::Bool=false
)
```

What would you like to spend time working on together?

Activity:

1. Use `supertypes(...)` and `print.(tt(...))` to inspect the type trees of ENA models and rotations
2. What is an abstract type and why do you think we see them in type trees so much?

`plot(model)` produces a plot with the following subplots:

- `(a)` an overall mean, which tells us the baseline everything compares against
- `(b)` and `(c)` rates of change for each connection across the x- and y-axes, which tells us what is *actually* being modeled by each axis. If you are coming up with the name of an axis after-the-fact, it's good to check your assumptions against these trends and make sure your name for the plot captures what is actually being modeled. If all the lines look grey, then it's hard to succinctly say *what* is being modeled on an axis: there's just a lot of noise
- Subsequent subplots show each subgroup on its own. It's good to compare these to the overall mean
- And the last subplots show how each pair of subgroups compare. Similar to the trend plots, these show you *what* is being modeled by the difference of the two groups. If everything looks grey, then it's hard to give a succinct description of that difference: there's just a lot of noise

Some differences from WebENA and rENA:

- Saturation shows *correlation strength* so we can tell at a glance when a model might be weak. In WebENA and rENA saturation is redundant (not necessarily bad) with line thickness, which shows magnitude of an effect
- Plots are mean centered by moving the mean of the plot, not by changing the underlying data. This preserves information that may or may not be useful for downstream analyses
- Plots are opinionated. Based on the model type, rotation type, and config, *I*, Mariah Knowles specifically, made the call for what the best way to plot that model is and made it the default. These defaults were chosen to be predictable and based only on the details you explicitly say yourself when you setup your model. This gives you the "right" plot without having to specify what "right" means each time
- A [known issue](https://github.com/snotskie/EpistemicNetworkAnalysis.jl/issues/11) is that the y-axis label can get cutoff when there are a lot of subplots

There are a lot of arguments for tweaking the default plot. For info on these, check the help text on `EpistemicNetworkAnalysis.plot` or [read the plot documentation](plots.md)

Some other useful functions, we can export an ENA model to a spreadsheet, we can save our plots to many image formats, and we can load in our own data from a CSV file:

```julia
to_xlsx("example.xlsx", model)

using Plots
p = plot(model)
savefig(p, "example.png")

using DataFrames
using CSV
data = DataFrame(CSV.File("example.csv"))
```

## Geometry; or, Picking the Right Rotation for the Job

Activity:

1. Put the image below on the board
2. Split the group in half: ethnographers and ethnographers-of-ethnographers
3. Ethnographers: How would you describe this image? Discuss and put sticky note labels right on the TV. Label everything you see
4. E-o-Es: Watch their discussion and take notes
5. When done, take the image away, leaving just the stickies
6. E-o-Es: What did you observe and what do you see going on up here?
7. Ethnographers: Watch their discussion and take notes
8. When done, Ethnographers: What did you observe?

![](map-capture.png)

All stories (or, the kind we tell) have a geometric intuition behind them:

- We talk about some set of events because they share some sense of similarity
- We move onto the next paragraph where we talk about a different set of events sharing a different similarity
- And we hang those paragraphs together into a larger narrative structure
- These narratives have a geometry to them: over here is this, and over here is this, and here's how the qualities change as we move around from point to point, from paragraph to paragraph

The ENA process as a whole has five steps:

1. Accumulate connection counts between qualitative codes over a sliding window. This gives a high dimensional matrix representation of one's discourse
2. Embed a network into that space as a way to approximately understand how the model weights change as one moves through the space
3. Reduce the dimensionality of that space to highlight one's features of interest. This is called dimension reduction, multidimensional scaling, or rigid body rotation in QE-land
4. Visualize that reduced space
5. Interpret the results

All these steps are agnostic to your sense of your narrative structure, except for the third, the dimension reduction step. It's in that step that the algorithm chooses which information to show and which to hide. Ideally, your dimension reduction lets you think on the plot the same way you think on the page&mdash;over here is this, and over here is this&mdash;, letting you (a) move through the story as feels natural while (b) retaining exactly the information you need to tell and test that kind of story

**Want a high-level view of your data?**

You need to see the most "stuff" of your story

You probably don't have a strong sense of a narrative structure just yet

Use `SVDRotation`, captures the most variance, captures the most "stuff," plot is the most likely to be legible at a glance

**Want to compare groups?**

You need to see what makes groups different from each other, and everything else is beside the point

If you have two groups, your story likely looks like a Venn diagram: Here's what group A has, here's what group B has, and here's what they have in common

Use a `MeansRotation` in that case, captures the most between-group variance, plot places the group means right on the x-axis

If you have more then two groups, then there's a few options.

If your groups have an innate order to them, like months in a year or phases of an intervention, then I'd like to sit down with you! I haven't found a convincing way to model these just yet

If your groups lack any innate order, like schools in a district, then your narrative structure might be hard to intuit. There's lots of ways to compare multiple groups! I recommend a "themes of difference" approach, where each paragraph is one theme that captures a lot of what makes multiple groups differ from one another

Use `MulticlassRotation` or `LDARotation`. `MulticlassRotation` retains more variance, so plots tend to be more legible and interpretable. `LDARotation` focuses on discrimination between groups, so it could be useful if you have downstream predictive analysis you want to run

**Want to relate your model to continuous variables?**

Like time, test scores, etc.

You need to see how the qualities of your discourse covary with your continuous variable, and everything else is beside the point

Your narrative structure probably moves through touch points in your continuous variable, such as low scores, medium scores, and high scores

If there is a well-defined and meaningful way to split your continuous variable into discrete groups, then consider telling your story by directly comparing those groups. If there aren't great ways to split your continuous variable into discrete groups, then consider using high dimensional clustering techniques to define groups of similar points, then tell the story that compares those groups

There are [many ways](https://link.springer.com/chapter/10.1007/978-3-031-31726-2_5) to think about stories that compare multiple groups. For some, a `FormulaRotation` is most appropriate:

- When you just want to see a general trend from low to high groups, without directly comparing groups
- When you don't have discrete groups and you just want to see a general trend of "high scores are like this" and "low scores are like this"
- When you need to control for confounding variables or nested data

`FormulaRotation` fits a linear regression to each dimension of the high dimensional space and uses that to find a trend line or "slope," which it runs the x-axis through. A good idea is to add a trajectory to the plot to help interpret the x-axis in terms of your continuous variable

**Have a theory about what makes your groups differ?**

Or any other statistical test really

You need to see the variance that is explained exactly in terms the codes/connections that make up your theory, not to *infer* those codes/connections from your plot

Your narrative structure likely looks like one of the narrative structures described above, but with an added paragraph of, "now let's test that idea directly."

Use `TopicRotation`, which let's you specify which codes/connections to force to the left and which to force to the right

**Need to test the generalizability or validation of your model?**

Here the goal isn't to tell a story, but to strengthen certain positivist evaluative criteria that apply to your work

Use `TrainedRotation` to create a new model that uses the same embedding and runs the same tests as an existing model

## Closeout

TODO