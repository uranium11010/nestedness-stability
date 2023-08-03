unitsize(2.5inch);

real[] Xs = {0, 0.19, 0.24, 0.29, 0.31, 0.39, 0.45, 0.51, 0.55, 0.58, 0.60, 0.66, 0.70, 0.75, 0.78, 0.83, 0.89, 0.93, 0.97, 1};
real[] Ys = {0, 0.03, 0.05, 0.07, 0.09, 0.15, 0.18, 0.20, 0.24, 0.26, 0.29, 0.33, 0.38, 0.45, 0.47, 0.53, 0.65, 0.71, 0.82, 1};
real[] extraXs = {0.06, 0.09, 0.16, 0.34, 0.47, 0.64, 0.72, 0.87};
real[] extraYs = {0.11, 0.13, 0.22, 0.33, 0.345, 0.36, 0.40, 0.425, 0.49, 0.545, 0.565, 0.59, 0.63, 0.68, 0.74, 0.755, 0.79, 0.85, 0.9, 0.92, 0.95, 0.985};
real x_low = 0.55, x_high = 0.58;
real y_low = 0.425, y_high = 0.45;
real x_ex = 0.92;
real y_ex = 0.17;
real d = 0.03;  // tick length
real opacity = 0.2;
pen grid_color = lightgray;
pen demons_color = royalblue;
pen ij_color = heavygreen;

assert(Xs.length == Ys.length);
int S = Xs.length - 1;
write("S = ", S);

real x = (x_low + x_high) / 2;
real y = (y_low + y_high) / 2;

// HELPERS

int get_x_idx(real x, real[] Xs) {
  int x_idx = 1;
  while (x_idx < S && Xs[x_idx] <= x) {
    ++x_idx;
  }
  return x_idx;
}

int get_y_idx(real y, real[] Ys) {
  int y_idx = S-1;
  while (y_idx > 0 && Ys[y_idx] >= y) {
    --y_idx;
  }
  return y_idx;
}

void draw_x_tick(string label, real x, real d, pen color=black) {
  draw(Label(label, BeginPoint), (x, 1+d) -- (x, 1), color);
}

void draw_y_tick(string label, real y, real d, pen color=black) {
  draw(Label(label, BeginPoint), (-d, y) -- (0, y), color);
}

// background
fill(unitsquare, white);

// grid lines
for (int i = 1; i < S; ++i) {
  draw((Xs[i], 0) -- (Xs[i], 1), grid_color);
  draw((0, Ys[i]) -- (1, Ys[i]), grid_color);
}
for (int i = 0; i < extraXs.length; ++i) {
  draw((extraXs[i], 0) -- (extraXs[i], 1), grid_color);
}
for (int i = 0; i < extraYs.length; ++i) {
  draw((0, extraYs[i]) -- (1, extraYs[i]), grid_color);
}

// linkage existence boundary
path boundary = (0, 0);

for (int i = 1; i <= S; ++i) {
  boundary = boundary -- (Xs[i], Ys[i-1]) -- (Xs[i], Ys[i]);
}
path no_link_region = boundary -- (1, 0) -- cycle;
fill(no_link_region);
// draw(boundary, grid_color);

// cell
draw((x_low, y_low) -- (x_low, y_high) -- (x_high, y_high) -- (x_high, y_low) -- cycle, ij_color);

// diagonal
draw((0, 0) -- (1, 1), demons_color);

// lower-left area
int x_idx = get_x_idx(x, Xs);
real Y = Ys[x_idx-1];
fill((x_low, 1) -- (x_low, Y) -- (x_high, Y) -- (x_high, 1) -- cycle, ij_color+opacity(opacity));
label("$j$", (x, 1), N, ij_color);
draw((Y, Y) -- (x_high, Y), demons_color);
path DLboundary = (Y, Y) -- (0, 0);
for (int i = 1; i < x_idx; ++i) {
  DLboundary = DLboundary -- (Xs[i], Ys[i-1]) -- (Xs[i], Ys[i]);
}
DLboundary = DLboundary -- cycle;
fill(DLboundary, demons_color+opacity(opacity));
int y_ex_idx = get_y_idx(y_ex, Ys);
real X_ex = Xs[y_ex_idx+1];
draw(Label("$w(y)$", Relative(0.45), align=Relative(E/2)), (y_ex, y_ex) -- (X_ex, y_ex), demons_color, Arrows);

// upper-right area
int y_idx = get_y_idx(y, Ys);
real X = Xs[y_idx+1];
fill((0, y_high) -- (X, y_high) -- (X, y_low) -- (0, y_low) -- cycle, ij_color+opacity(opacity));
label("$i$", (0, y), W, ij_color);
draw((X, X) -- (X, y_low), demons_color);
path DLboundary = (X, X) -- (1, 1);
for (int i = S-1; i > y_idx; --i) {
  DLboundary = DLboundary -- (Xs[i+1], Ys[i]) -- (Xs[i], Ys[i]);
}
DLboundary = DLboundary -- cycle;
fill(DLboundary, demons_color+opacity(opacity));
int x_ex_idx = get_x_idx(x_ex, Xs);
real Y_ex = Ys[x_ex_idx-1];
draw(Label("$h(x)$", Relative(0.6), align=Relative(E/2)), (x_ex, x_ex) -- (x_ex, Y_ex), demons_color, Arrows);

// left ticks
draw_y_tick("0", 1, d, demons_color);
draw_y_tick("$y_\mathrm{min}$", Y, d, demons_color);
draw_y_tick("$y$", y_ex, d, demons_color);
draw_y_tick("1", 0, d, demons_color);
draw((0, Y) -- (Y, Y), dashed+demons_color);
draw((0, y_ex) -- (y_ex, y_ex), dashed+demons_color);

// top ticks
draw_x_tick("0", 0, d, demons_color);
draw_x_tick("$x_\mathrm{min}$", X, d, demons_color);
draw_x_tick("$x$", x_ex, d, demons_color);
draw_x_tick("1", 1, d, demons_color);
draw((X, 1) -- (X, X), dashed+demons_color);
draw((x_ex, 1) -- (x_ex, x_ex), dashed+demons_color);

// point
label("$(i, j)$", (x, y), N+E, ij_color);

// border
draw(unitsquare);

// background grid (for guiding drawing)
// for (int i = 1; i < 10; ++i) {
//   draw((0, i/10) -- (1, i/10), green);
//   draw((i/10, 0) -- (i/10, 1), green);
// }