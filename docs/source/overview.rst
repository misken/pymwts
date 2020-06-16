===============================================
Project overview
===============================================


Service systems such as hospitals and call centers face tough
staff scheduling problems. Staffing needs often vary significantly
by time of day and day of week. The use of flexible scheduling
practices such as the use of different shift lengths, start times,
and the number of hours worked per week, while making it easier to
match variable staffing requirements, can lead to large and complex staff
scheduling problems. 

There are a bunch of different flavors of staff scheduling problems.
For example, the ongoing task of creating four week schedules for a pool of nurses working in a set of nursing units is quite different than trying to determine which mix of shift lengths provides the best match to hourly staffing needs in a post-anesthesia care unit. We call this
second type of problem a *tactical scheduling analysis* problem and it
is these types of problems for which this model can play an important role.

Some history
------------

An early version of this model was used in numerous real tactical
scheduling problems back when one
of the developers was a practicing operations analyst for a few large
healthcare systems. After leaving industry for academia, the model and
description of its use
was published as a journal article. 

Isken, Mark W. "An implicit tour scheduling model with applications in healthcare." *Annals of Operations Research* 128.1-4 (2004): 91-109.

The original model is a one-week tour scheduling model (tour scheduling
is a combination of shift and days worked scheduling) and while useful,
has several shortcomings related to its single week nature::

- it is difficult to model tours in which the number of days worked varied over the weeks,
- it is difficult to model different weekend policies,
- it is difficult to construct sample multiweek schedules from the one week solution.

We use the word "difficult" since we would do all kinds of gymnastics to try to
create one week schedules that we could then manually postprocess and combine to
create example multiweek schedules for our projects.

Even though these kinds of models are not appropriate for ongoing operations schedule creation and
management, it is highly desirable to be able to produce multiweek schedules to 
show decision makers that realistic schedules can be constructed and visualized for various scheduling policies.

Problem size challenges
-----------------------

Why was the original model a one-week model? Well, these are combinatorial
optimization problems and multiweek tour scheduling problems get 
huge very quickly as one considers increased scheduling flexibility. In
traditional tour scheduling models, there is one variable per allowable tour.
The number of variables can grow into the millions or billions (or more) for quite
realistic scenarios for scheduling flexibility. Even for one
week models, the problems can get quite large when modeled in the
traditional way. Our one-week model took a different approach known as
*implicit modeling* to get around the problem size explosion. Again,
see the paper above for all the details.

Extending the implicit model to handle multiple weeks
------------------------------------------------------

So, the original model sat around for a bunch of years, but I'd always
planned on extending it to multiple weeks. Now that a usable multiweek
model is complete and a journal article is under review (SOON), the model code is being released as an open source
project called `pymwts`. You can see a `preprint version of the paper here (COMING SOON) <>`_.

About the code
--------------

The model was implemented in Python using the open source optimization
modeling language `Pyomo <http://www.pyomo.org/>`_. It can be used
with any solver with mixed integer programming capabilities such as
`CBC <>`_, `GLPK <>`_ or `Gurobi <>`_.

This is what I'd call "research quality" code. It's pretty well commented
and the formulation as presented in our paper intentionally matches
the way the model was implemented in Pyomo. However, I'm sure there's
much room for improvement, especially as the code evolved over the
years as I worked on it here and there.

It's got a simple command line interface.

About the sample data files
---------------------------

If you have some experience with using algebraic modeling languages you
know that one of the advantages is the ability to segregate optimization
model logic from the instance specific data. In Pyomo, our model
is implemented as an `AbtractModel` and the data files are written in
standard `AMPL compatibile DAT format <https://ampl.com/BOOK/CHAPTERS/24-refman.pdf>`_. Pyomo has `several other options for
specifying data <https://pyomo.readthedocs.io/en/stable/working_abstractmodels/instantiating_models.html>`_ for an optimization model instance but DAT files are a simple, portable,
transparent way for us to provide example problem instances. 

The `pymwts` project contains some code for generating DAT files from
YAML formatted input files but more work is needed to make these
generally usable.

