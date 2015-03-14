# Introduction #

This wiki is intended to help users get started with the code. If you are interested in the code and want to keep up to date, please join the [google group](https://groups.google.com/group/lily-pad). If you think there is an additional step worth mentioning, please add it. If you need further help, please contact the admin. If you find issues with the code, please note them in the issue tracker.

For documentation on the algorithms and project which have used Lily Pad, please see the LilyPadDocumentation.


# Getting Started #

The code is written in Processing and uses Git for version control. Before working with Lily Pad you will need to get these two programs (both of these are cross platform, so get the flavor you need) :
  * Get a copy of [Processing](http://processing.org/). They have an excellent [wiki](http://wiki.processing.org/).
  * Get a copy of [Git](http://git-scm.com/). There are many nice [introductions](http://git-scm.com/documentation).

Now you can download the Lily Pad program.  Using the google code website to make a clone is **not recommended**. It is pretty confusing. Instead, at the command line type:

`$ git clone https://code.google.com/p/lily-pad/ LilyPad`

That's it! You're ready to run.


# Running #

Open the code in Processing and click the _play_ button. The default set-up (a cylinder in a free stream) will run automatically. You can drag the cylinder around in the default set-up as well.

The main Processing file is `LilyPad.pde` and this file controls what happens when you click play. Each of the files in the `LilyPad` folder has an example code section commented out at the top. By replacing the contents of the `LilyPad.pde` file with these example code sections you can run a wide variety of simulations and explore the types of geometries, flow conditions, and visualizations available.


# Contributing #

However, soon enough you will probably want to do something which isn't coded up yet. Great! This is an active and growing project and we love seeing how people are using it.

First, go to the [google group](https://groups.google.com/group/lily-pad), log in with your google ID, and click on `Join Group` to be added to the project group. This will keep you up to date with code changes. Also, you need to set your name and email address within the git environment. A good list of first time set set-up stuff is [here](http://git-scm.com/book/en/Getting-Started-First-Time-Git-Setup).

And as you contribute, don't worry too much about making a mistake or putting up ugly code. We can always fix/clean it. Even better, we can suggest how you can fix/clean it yourself. IF YOU HAVE PROBLEMS, go to the [google group](https://groups.google.com/group/lily-pad) and ask a question!


## Best practises ##

The basic guideline is to add code which is general and friendly enough to be useful for other users. This means:
  * Creating and pushing test case files (such as `AudreyTest.pde`) is especially useful since it gives everyone a chance to run and play with a new test case out-of-the-box!
  * Let users run your code through **high-level** interfaces, and only dig down into the guts of the methods if they choose.
  * Reuse as much existing code as possible. i.e. add a new capability by making a new class which **extends** a current class.
  * If you change low level routines, do so in a way which does not break compatibility. i.e. trigger your change with an **optional** argument.
  * Completely document and test all changes before pushing up.
  * Add a minimum working example at the top of any new files.
  * Use the [Code Review](http://code.google.com/p/support/wiki/CodeReviews) feature to request (or assign) a review of the new code you have committed.


## Pushing changes up ##

Once you have made (and tested and `commit`ted; see the introduction to git pages above) a set of changes you want to share with others, you need to `push` them up to the repository. To do this you simple type

`git push`

Supply your gmail username and your personal google code password, found [here](http://code.google.com/hosting/settings), and you're done! Note that some git clients don't allow copy-paste or give you any indication that they are accepting a password (as in, they don't even fill up with `*` characters). Dont worry - it hasn't hung up - just type carefully and hit enter.

Good luck!