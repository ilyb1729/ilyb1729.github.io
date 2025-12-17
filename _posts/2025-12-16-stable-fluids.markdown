---
layout: post
title: Stable Fluids
date: 2025-12-16 00:00:00 -0700
categories: [] # Optional categories
tags: [] # Optional tags
---

I have recently been working on an implementation of the stable fluids algorithm and thought it might be good to collect some of the process and my thoughts since its been about a year since I last attempted to implement a fluid sim for CS488 (Computer Graphics).

# Intuition

## Definitions/Conventions

For this, I will stay in 2D to keep notation simple and because my sim is implemented for 2D only (it's easy to generalize to 3D).

We represent velocity as the vector $\mathbf{u} = (u, v)$, time step as $\Delta t$, viscosity as $\nu$, density as $\rho$, and external force as $\mathbf{f}$.

In Stam, they present Navier-Stokes as
\\[
\nabla \cdot \mathbf{u} = 0
\\]
and
\\[
\frac{\partial \mathbf u}{\partial t} = -(\mathbf{u} \cdot \nabla)\mathbf{u} - \frac{1}{\rho} \nabla p + \nu \nabla^2\mathbf{u} + \mathbf{f}.
\\]
They introduce an operator $\mathbf{P}$ that lets you project a vector field into its divergence free part based on the Helmholtz-Hodge decomposition (not sure what this is). Then, we get
\\[
\frac{\partial \mathbf u}{\partial t} = \mathbf{P}\left(-(\mathbf{u} \cdot \nabla)\mathbf{u} + \nu \nabla^2\mathbf{u} + \mathbf{f}\right).
\\]
In my code, I am using something called a MAC (marker-and-cells) grid. What this means is that vertical velocities are stored on the center horizontal edges and horizontal velocities are stored on the center of vertical edges. This supposedly helps with artifacts from axis aligned velocities and will change what the discretized formulas look like. For a grid of size $m \times n$, this means that the grid storing $u$ velocities is $m+1 \times n$ and the grid storing $v$ velocities is $m \times n + 1$.

## Steps

Now, we introduce a high level view of what each of the steps does. I referenced my approach from the stable fluids paper and during my design phase. I wanted to read a textbook on fluid sims to get my information from a more reliable source, but my main goal during this project was not the accuracy of simulation so I did what makes sense swiftly.

### Applying forces

Consider the $\mathbf{f}$ on the RHS of Navier-Stokes. This is pretty self explanatory, you apply the external forces otherwise there will be no velocity.

### Advection

Now we consider the $-(\mathbf{u} \cdot \nabla) \mathbf{u}$ part. It was confusing for me at first since I thought dotting a field and an operator doesn't type check, but if you just accept the abuse of notation, we can see that
\\[
u \cdot \nabla \equiv u_1 \frac{\partial u}{\partial x} + u_2 \frac{\partial u}{\partial y} \equiv D_u.
\\]
So this is telling us that velocity moves along its own field (direction derivative of velocity field along the velocity field).

Stam shows how to compute the transfer of velocity through backtracking a simulated particle through the velocity field. They choose to backtrack through time instead of forward through time because it follows the implicit approach which supposedly makes it more numerically stable (implicit means something like the result is based on the outcome?).

Anyways, the idea behind the backtracking is that if our particle has position $\mathbf{x}$, then we want to align the position of the particle with the diffeq
\\[
\frac{\partial \mathbf{x}}{\partial t} = \mathbf{u}(x, t),
\\]
backtrack through time, find the velocity at that point in space/time, and set the velocity at the particle's initial starting time to the newly computed velocity.

#### Aside: Runge-Kutta

Runge-Kutta allows you to solve for diffeqs of the form above. RK2 is the second order version and we present one possible form here. The convergence conditions seem pretty complicated, so I just ripped what works already.

What we do is compute the velocity at the end point and move back based on the velocity
\\[
x_{middle} = x_{end} - \frac{\Delta t}{2} v(x_{end}).
\\]
Then we treat this velocity $v(x_{middle})$ as the estimated velocity over the interval from $0$ back in time to $-\Delta t$ and estimate
\\[
x_{beginning} = x_{end} - \Delta t \cdot v(x_{middle}).
\\]
This gives us the initial velocity $v(x_{beginning})$.

### Diffusion

Next we have $\nu \nabla^2 \mathbf{u}$. This just lets velocity spread from/to its neighbors depending on the viscosity of the fluid. Apparently it uses the second order partials since viscous force is related to the divergence of shear and shear is a first order derivative. I don't understand the physics as I live in a computer, not in the real world.

If we use a finite difference approximation of a derivative, we can approximate
\\[
\nabla^2 u = \frac{u_{up} - 2 u_{center} + u_{down} + u_{right} - 2 u_{center} + u_{left}}{\ell^2}.
\\]
$\ell$ is the length of a cell, or $1$ in my simulation to simplify things.

Since we have
\\[
\frac{\partial u}{\partial t} = \nu \nabla^2 \mathbf{u},
\\]
we can rewrite in the implicit form to get
\\[
\frac{\mathbf u^{n+1} - \mathbf u^n}{\Delta t} = \nu \nabla^2 \mathbf{u^{n+1}} \implies \mathbf u^{n+1} - \Delta t \cdot \nu \nabla^2 \mathbf u^{n+1} = \mathbf u^n .
\\]
This gives us a system of equations where the unknown are the velocities in $\mathbf u^{n+1}$.

#### Aside: Gauss-Siedel

Give a system in the form $Ax = b$, we can rewrite $A = L + U$ where $L$ is a lower triangular matrix and $U$ is a strictly upper triangular matrix (is this the right term? no nonzero elements on the diagonal). Then, we can rewrite
\\[
(L + U)x = b \implies Lx = b - Ux.
\\]
As an approximation, we can let the LHS be future $x^{n+1}$ and RHS be the previous $x$, so we get the iteration
\\[
Lx^{n+1} = b- Ux^{n}.
\\]
We can now solve this with standard triangular matrix techniques. Convergence is not obvious. maybe something to learn in the future.

### Projection

Handwavily, we start by considering the $-\frac{1}{\rho} \nabla p$ part of Navier-Stokes. Discretizing over our time step, we get $-\frac{\Delta t}{\rho}\nabla p$. Since we are assuming constant density for our simulation, we can botch the physical meaning and absorb the constants into $\nabla p$ due to linearity of partial derivatives ($p$ I think means pressure?) to get $\nabla \tilde p$ and we get
\\[
u^{n+1} = u^* - \nabla \tilde{p}.
\\]
We must have that $\nabla u = 0$ from Navier-Stokes for our new value $u^{n+1}$, so
\\[
\nabla (u^* - \nabla \tilde{p}) = 0 \implies \nabla u^* = \nabla^2 \tilde{p}.
\\]
We know $u^*$ from our previous steps and $\tilde p$ is unknown. For this system to make sense, we evaluate it for each cell in our MAC grid. The gradient can be discretized as
\\[
\nabla u = \frac{u_{up} - u_{down} + u_{right} - u_{left}}{\ell}
\\]
and the laplacian can be discretized with a derivative of a derivative as

$$(\nabla^2 p)_{i, j} = \frac{p_{i+1, j} - 2p_{i, j} + p_{i-1, j} + p_{i, j+1} - 2p_{i,j} + p_{i,j-1}}{\ell^2}$$

This gives us a system of equations with a variable per cell. We can solve this using the same Gauss-Siedel approach from the last step.

### Boundary Conditions

One thing that is left out of all of my math is that in the code you have to deal with boundaries that are really pesky. Actually, I'm fairly certain that I got cooked by my boundary conditions. Anyways, the two most common seem to be slip and non-slip boundary conditions.

In slip, you make all normal velocities zero. This means zeroing out the velocities on the boundary. So in a MAC grid, all of the horizontal velocities on vertical borders and vertical velocities on horizontal borders get zeroed out.

In non-slip, you make all normal and tangential velocities zero. This means doing the same thing as slip and also for vertical borders, you make the vertical velocities outside the border the negative of vertical velocities mirrored across the border and vice versa for horizontal borders. This ensures that when you perform bilinear interpolation on any point on the border, you will always get 0 velocity in both $u$ and $v$ directions.

On a related note, in systems like the one for projection where it relies on velocities that might be zeroed out or not in the grid, I just left them out by cutting out a buffer of one cell from the system that I solve. I saw some similar approaches in other people's implementations.

# Software Approach (Testing?)

There were quite a few things that tripped me up last time I attempted a fluid sim:

1. Off by one errors with MAC staggered grids
2. Slow and poor visualization of what was happening
3. Starting in 3D
4. Trying to test too many steps at the same time

Basically all these boil down to my strategy for testing and debugging. Recently I have been buying more into my friends idea that your ability to iterate and test your mental model against the world are essential. In fact, I think it might be the greatest indicator of your ability to code.

The first major change I made was implementing this in 2D. Not only did this make the code simpler so that I could focus on the process, it also made the visualization much easier. I made sure that the first thing I set up were helper functions that let me visualize the $u$, $v$ velocities and the speed and quickly figured out how to generate velocity fields. This way, the rest of my testing was much easier and I could get instant feedback on what is going wrong by setting up different situations and changing the type of visualization.

I also put in enough effort to write unit tests for my staggered grid implementation! The trauma was so bad last time I tried to implement one with one off errors, but I noticed that my mental model was much better and the tests helped build confidence in my implementation.

For the rest of my tests, I thought the most reasonable way was to implement one piece of the simulation at a time, run some tests, and validate that the fluid is behaving in a reasonable way. This way, I am able to tighten my compile, test, revise loop and keep my mental model in sync with what the code is doing. I found that this worked well enough for my goal of having something workable that I can optimize for this project.

As for possible further testing, I was considering snapshot testing, but it felt suspicious to do so since basically any code change will either cause me to need to update snapshot (no matter what values I am saving) or it would be hard to separate out a breaking change and a proper change by setting a threshold epsilon diff since we are dealing with floating points.

I am still a noob at theorem prover languages, but I contemplated if these languages would help with writing a fluid sim. The issue is that since the algorithms are just an approximation of a physical phenomena, the specification of the algorithm and the actual implementation would be almost the same thing. Perhaps it might be possible to prove some properties about the approximation properly following Navier-Stokes or the numerical stability of the algorithm. Both of these seems incredibly hard though. I could be completely wrong about everything here as I have not written much code in this area.

# Reflection on Coding

I started to work on the abstractions after figuring out what the general high level ideas to the level of diffusion spreads particles, projection enforces incompressibility, etc. This was a good step in the direction of designing more carefully before coding, but there were still a few times where I had to rethink how I was storing my data and what my abstractions for grid vs simulation ended. This is partially since I am inexperienced in sims still; I had to do much less backtracking compared to my first attempt and I anticipate that next time I write a fluid sim it will only get much easier. However, I did not understand the computation in enough detail before implementation. I should have thought about how the data would be transformed. For example, I guessed on how to store the velocities: two separate 2D grids with one velocity component or a singular 2D array with two components per cell. I got lucky that the access patterns computing one direction at a time and also the different sizes required to implement a MAC grid ($(n+1) \times n$ vs $n \times (n+1)$) made the first choice better, but next time I should be able to figure this out before having to guess and write code.

Additionally, I tried to not use GPT to write my code. This is because coding is something I am trying extremely hard to work on and I spend enough time at work practicing on my AI workflows so I think this is the best way to use my time. I did slip up a little bit when I got frustrated and tried to brute force some progress with GPT, but I woke up the next day and fixed the slop that GPT produced. But there is still progress to be made to completely separate myself from GPT (and determining if this is efficient to my development for the future).

One cost of not using GPT is that documenting my work after doing all of the math and coding in my head is extremely painful. Documenting would have made my code easier and doing documentation in parallel would make the documentation easier. Will do that next time. I am used to getting away with it since AI helps a lot and I don't write as mathematical stuff normally.

As a side note, still struggling with rust a little bit, but things are slowly starting to make sense. Since Gauss-Siedel works best in place and copying large arrays is prohibitively expensive, I was forced to use mutable references more than I was used to so I got to understand shared reference restrictions better. I prefer to learn things from the fundamentals but I heard the borrow checker is extremely complicated. If you have any suggestion for a lower level explanation for how the borrow checker works, please let me know!

# Future Work

I believe that the current results look somewhat reasonable (I have to figure out how to embed gifs into this :/) as it matches what I expect the equations to solve for. Though, there are some things that look wrong from a physical standpoint of when I observe water. This will be something I may circle back to though, since the main purpose of this project was to have a small piece of nontrivial code so that I can try out some CPU optimization.

I will maybe have a followup post about some things I do as I try to get my simulator to run above 2 FPS!

Feel free to email me rcheng1729@gmail.com if you have any ideas or comments.
