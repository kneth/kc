/*
  Golden-output test for one kc model in one mode.

  usage: runner KC MODEL MODE EXPECTED_DIR

  Runs "KC -mMODE < MODEL" in a fresh directory named MODEL.mMODE and
  compares every generated file with EXPECTED_DIR/MODEL.mMODE/. In mode 3
  the generated model is also compiled with kksolver and simulated when the
  expected output contains sim.last (the final state of the simulation).

  Numbers are compared with a relative tolerance and timestamps are ignored,
  so results are portable across compilers and platforms. Entries listed in
  KC_SKIP_LIST are skipped: whole cases ("model.mN") or single files
  ("model.mN/file").

  Set KC_UPDATE=1 to (re)write the expected output instead of comparing.
*/

#include <ctype.h>
#include <dirent.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

#define SKIP 77        /* CTest SKIP_RETURN_CODE */
#define REL_TOL 1e-6
#define SIM_REL_TOL 1e-4 /* step size control differs slightly between compilers */
#define ABS_TOL 1e-10

static bool update;

/* snprintf that treats truncation as a fatal error. */
static void format(char *buf, size_t size, const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  int n = vsnprintf(buf, size, fmt, ap);
  va_end(ap);
  if (n < 0 || (size_t)n >= size) {
    fprintf(stderr, "runner: path or command too long\n");
    exit(2);
  }
}

static char *slurp(const char *path) {
  FILE *f = fopen(path, "rb");
  if (!f) return NULL;
  fseek(f, 0, SEEK_END);
  long n = ftell(f);
  rewind(f);
  char *s = malloc(n + 1);
  s[fread(s, 1, n, f)] = '\0';
  fclose(f);
  return s;
}

/* Matches a ctime(3) timestamp such as "Fri Oct 30 14:48:31 2026". */
static bool is_date(const char *s) {
  const char *pat = "AAA AAA _9 99:99:99 9999";
  for (; *pat; pat++, s++) {
    if (*pat == 'A' && !isalpha((unsigned char)*s)) return false;
    if (*pat == '9' && !isdigit((unsigned char)*s)) return false;
    if (*pat == '_' && *s != ' ' && !isdigit((unsigned char)*s)) return false;
    if (strchr(" :", *pat) && *s != *pat) return false;
  }
  return true;
}

static void mask_dates(char *s) {
  for (; *s; s++)
    if (is_date(s)) memset(s, '#', 24);
}

static bool starts_number(const char *s, const char *begin) {
  if (s > begin && (isalnum((unsigned char)s[-1]) || s[-1] == '_' || s[-1] == '.'))
    return false;
  return isdigit((unsigned char)*s) || strchr("+-.", *s) || !strncmp(s, "nan", 3);
}

/* Compares text, treating numbers as equal within tolerance. Returns the
   1-based line of the first difference, or 0 when equal. */
static int compare(const char *a, const char *b, double rel_tol) {
  const char *a0 = a, *b0 = b;
  int line = 1;
  while (*a || *b) {
    char *ea, *eb;
    if (starts_number(a, a0) && starts_number(b, b0)) {
      double x = strtod(a, &ea), y = strtod(b, &eb);
      if (ea != a && eb != b) {
        bool same = (isnan(x) && isnan(y)) || x == y ||
                    fabs(x - y) <= ABS_TOL + rel_tol * fmax(fabs(x), fabs(y));
        if (!same) return line;
        a = ea;
        b = eb;
        continue;
      }
    }
    if (*a != *b) return line;
    if (*a == '\n') line++;
    a++;
    b++;
  }
  return 0;
}

static bool skipped(const char *entry) {
  char *list = slurp(KC_SKIP_LIST), *tok;
  bool found = false;
  if (!list) return false;
  for (tok = strtok(list, "\n"); tok && !found; tok = strtok(NULL, "\n"))
    found = strcmp(tok, entry) == 0;
  free(list);
  return found;
}

static int run(const char *cmd) {
  int status = system(cmd);
  return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
}

/* Compiles the generated model with kksolver, simulates it and stores the
   last line of the simulation output in DIR/sim.last. */
static bool simulate(const char *dir) {
  char cmd[4096];
  format(cmd, sizeof cmd,
           "mkdir -p %s/sim && cd %s/sim && "
           "'%s' %s -I.. -I'%s' '%s/kksolver.c' '%s' -lm -o sim 2> cc.log && "
           "./sim > /dev/null 2> sim.log && tail -n 1 kinwrkda > ../sim.last",
           dir, dir, KC_CC, KC_SIM_CFLAGS, KC_SOLVERS, KC_SOLVERS, KC_SOLVERS_LIB);
  return run(cmd) == 0;
}

static int check_file(const char *dir, const char *expdir, const char *name) {
  char path[1024], exp_path[1024];
  format(path, sizeof path, "%s/%s", dir, name);
  if (skipped(path)) return 0;
  format(exp_path, sizeof exp_path, "%s/%s", expdir, name);

  char *got = slurp(path), *want = slurp(exp_path);
  int rc = 0;
  mask_dates(got);
  if (update) {
    FILE *f = fopen(exp_path, "w");
    fputs(got, f);
    fclose(f);
  } else if (!want) {
    fprintf(stderr, "%s: no expected output\n", path);
    rc = 1;
  } else {
    mask_dates(want);
    int line = compare(got, want, strcmp(name, "sim.last") ? REL_TOL : SIM_REL_TOL);
    if (line) {
      fprintf(stderr, "%s: differs from expected at line %d\n", path, line);
      rc = 1;
    }
  }
  free(got);
  free(want);
  return rc;
}

/* Every expected file must have been generated (sim.last is checked apart). */
static int check_missing(const char *dir, const char *expdir) {
  DIR *d = opendir(expdir);
  struct dirent *e;
  char path[1024];
  int failures = 0;
  while (d && (e = readdir(d))) {
    if (e->d_name[0] == '.' || !strcmp(e->d_name, "sim.last")) continue;
    format(path, sizeof path, "%s/%s", dir, e->d_name);
    if (access(path, R_OK) != 0) {
      fprintf(stderr, "%s: expected file was not generated\n", path);
      failures++;
    }
  }
  if (d) closedir(d);
  return failures;
}

int main(int argc, char *argv[]) {
  if (argc != 5) {
    fprintf(stderr, "usage: runner KC MODEL MODE EXPECTED_DIR\n");
    return 2;
  }
  const char *kc = argv[1], *model = argv[2], *mode = argv[3];
  const char *base = strrchr(model, '/') ? strrchr(model, '/') + 1 : model;
  char dir[512], expdir[1024], cmd[4096], sim_last[1024];
  update = getenv("KC_UPDATE") != NULL;

  format(dir, sizeof dir, "%s.m%s", base, mode);
  if (skipped(dir)) return SKIP;
  format(expdir, sizeof expdir, "%s/%s", argv[4], dir);
  format(sim_last, sizeof sim_last, "%s/sim.last", expdir);
  if (update) mkdir(expdir, 0755);

  format(cmd, sizeof cmd, "rm -rf %s && mkdir %s && cd %s && '%s' -m%s < '%s' > stdout 2> stderr",
           dir, dir, dir, kc, mode, model);
  int status = run(cmd);
  if (status != 0) {
    fprintf(stderr, "%s: kc failed with status %d\n", dir, status);
    return 1;
  }

  bool want_sim = !strcmp(mode, "3") && (update || access(sim_last, R_OK) == 0);
  bool sim_ok = want_sim && simulate(dir);
  if (want_sim && !sim_ok && !update) {
    fprintf(stderr, "%s: simulation failed\n", dir);
    return 1;
  }

  int failures = check_missing(dir, expdir);
  DIR *d = opendir(dir);
  struct dirent *e;
  while ((e = readdir(d))) {
    char path[1024];
    struct stat st;
    format(path, sizeof path, "%s/%s", dir, e->d_name);
    if (stat(path, &st) == 0 && S_ISREG(st.st_mode) && (sim_ok || strcmp(e->d_name, "sim.last")))
      failures += check_file(dir, expdir, e->d_name);
  }
  closedir(d);
  return failures ? 1 : 0;
}
