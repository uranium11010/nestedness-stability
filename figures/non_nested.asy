size(2.5inch);

int R = 10;
int C = 6;
real cell_w = 1;
real cell_h = 1;
real figintvl = 1.4 * C * cell_w;

path cell = scale(cell_w, cell_h) * unitsquare;
path box = scale(C, R) * cell;

void draw_nested(int[] Xs, int[] Ys, transform t) {
  filldraw(t * box, white);
  assert(Xs.length == Ys.length);
  assert(Xs[0] == 0 && Xs[Xs.length-1] == C);
  assert(Ys[0] == 0 && Ys[Ys.length-1] == R);
  path boundary = (0, 0);
  for (int i = 1; i < Xs.length; ++i) {
    boundary = boundary -- (Xs[i] * cell_w, Ys[i-1] * cell_h) -- (Xs[i] * cell_w, Ys[i] * cell_h);
  }
  path no_link_region = boundary -- (C * cell_w, 0) -- cycle;
  fill(t * no_link_region);
}

void draw_nested_boundary(int[] Xs, int[] Ys, transform t) {
  path boundary = (0, 0);
  for (int i = 1; i < Xs.length; ++i) {
    boundary = boundary -- (Xs[i] * cell_w, Ys[i-1] * cell_h) -- (Xs[i] * cell_w, Ys[i] * cell_h);
  }
  draw(t * boundary, red+dashed);
}

struct Movement {
  pair begin, end;
  void operator init(pair begin, pair end) {
    this.begin = begin;
    this.end = end;
  }
}

void draw_move_cell(Movement movement, transform t) {
  filldraw(t * shift(scale(cell_w, cell_h) * movement.begin) * cell, drawpen=gray);
  filldraw(t * shift(scale(cell_w, cell_h) * movement.end) * cell, white, gray);
}

void draw_movement_arrow(Movement movement, transform t) {
  draw(t * scale(cell_w, cell_h) * (shift(0.5, 0.5) * movement.begin -- shift(0.5, 0.5) * movement.end), gray, Arrow);
}

int[] Xp = {0, 1, 2, 3, 4, 5, 6};
int[] Yp = {0, 1, 2, 4, 5, 7, 10};
Movement[] movements = {Movement((1, 2), (2, 1)), Movement((0, 1), (1, 0)), Movement((3, 5), (4, 3)), Movement((4, 7), (5, 4))};
draw_nested(Xp, Yp, identity);
for (Movement movement : movements) {
  draw_move_cell(movement, identity);
}
draw_nested_boundary(Xp, Yp, identity);
for (Movement movement : movements) {
  draw_movement_arrow(movement, identity);
}

void draw_deleted_cell(pair deleted, transform t) {
  draw(t * shift(scale(cell_w, cell_h) * deleted) * cell);
  draw(t * scale(cell_w, cell_h) * (deleted -- shift(1, 1) * deleted));
  draw(t * scale(cell_w, cell_h) * (shift(1, 0) * deleted -- shift(0, 1) * deleted));
}

int[] Xtilde = {0, 2, 3, 5, 6};
int[] Ytilde = {0, 1, 3, 4, 10};
pair[] deleted_cells = {(0, 1), (1, 2), (3, 3), (3, 5), (4, 4), (4, 7), (5, 5), (5, 6)};
draw_nested(Xtilde, Ytilde, shift(figintvl, 0));
for (pair deleted : deleted_cells) {
  draw_deleted_cell(deleted, shift(figintvl, 0));
}
draw_nested_boundary(Xtilde, Ytilde, shift(figintvl, 0));

label("$E' \to E$", (0.5 * C * cell_w, R * cell_h), N);
label("$\tilde E \to E$", (figintvl + 0.5 * C * cell_w, R * cell_h), N);