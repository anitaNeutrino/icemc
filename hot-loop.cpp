#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dlfcn.h>
#include <stdlib.h>
#include <string.h>

# include "TCanvas.h"
# include "TLine.h"
# include "TApplication.h"
# include "TSystem.h"

#include "game.h"
#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
#include <time.h>

#define NK_INCLUDE_FIXED_TYPES
#define NK_INCLUDE_STANDARD_IO
#define NK_INCLUDE_STANDARD_VARARGS
#define NK_INCLUDE_DEFAULT_ALLOCATOR
#define NK_IMPLEMENTATION
#define NK_XLIB_IMPLEMENTATION
#include "nuklear.h"
#include "nuklear_xlib.h"
#define DTIME           20
#define WINDOW_WIDTH    800
#define WINDOW_HEIGHT   600

#define UNUSED(a) (void)a
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) < (b) ? (b) : (a))
#define LEN(a) (sizeof(a)/sizeof(a)[0])


typedef struct XWindow XWindow;
struct XWindow {
    Display *dpy;
    Window root;
    Visual *vis;
    Colormap cmap;
    XWindowAttributes attr;
    XSetWindowAttributes swa;
    Window win;
    int screen;
    XFont *font;
    unsigned int width;
    unsigned int height;
    Atom wm_delete_window;
};

// struct nk_context *gctx = &xlib.ctx;
struct xlibstruct *gxlib = &xlib;

static void
die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputs("\n", stderr);
    exit(EXIT_FAILURE);
}

static void*
xcalloc(size_t siz, size_t n)
{
    void *ptr = calloc(siz, n);
    if (!ptr) die("Out of memory\n");
    return ptr;
}

static long
timestamp(void)
{
    struct timeval tv;
    if (gettimeofday(&tv, NULL) < 0) return 0;
    return (long)((long)tv.tv_sec * 1000 + (long)tv.tv_usec/1000);
}

static void
sleep_for(long t)
{
    struct timespec req;
    const time_t sec = (int)(t/1000);
    const long ms = t - (sec * 1000);
    req.tv_sec = sec;
    req.tv_nsec = ms * 1000000L;
    while(-1 == nanosleep(&req, &req));
}


char SO_LOCATION[256] = "";

char *errstr;

struct game {
    void *handle;
    ino_t id;
    struct game_api api;
    void *state;
};

static void game_load(struct game *game, bool bInteractive)
{
    struct stat attr;
    if (stat(SO_LOCATION, &attr) == 0) {
      printf("%s exist\n", SO_LOCATION);
      if (game->handle) {
        // usleep(500000000);
        game->api.unload(game->state);
        for (int i = 0; i < 1; i++){
          dlerror(); // Clear the state.
          int eret = dlclose(game->handle);
          if (eret) {
            printf("dlclose eret: %d\n", eret);
            printf ("dlclose err message: (%s)\n", dlerror());
          } else printf("dlclose successful: %d\n", eret);
        }
      }
      // usleep(200000);

      dlerror();
      void *handle = dlopen(SO_LOCATION, RTLD_NOW | RTLD_LOCAL);
      printf("handle: %p for %s\n", handle, SO_LOCATION);

      if (handle) {
        game->handle = handle;
        game->id = attr.st_ino;
        const struct game_api *api = (game_api *) dlsym(game->handle, "GAME_API");
        if (api != NULL) {
          game->api = *api;
          if (game->state == NULL) game->state = game->api.init(bInteractive);
          game->api.reload(game->state);
        } else {
          printf("api == NULL\n");
          dlclose(game->handle);
          game->handle = NULL;
          game->id = 0;
        }
      } else {
        printf("Couldn't get a handle\n");
        game->handle = NULL;
        game->id = 0;
      }
    } else {
      printf("%s doesn't exist\n", SO_LOCATION);
    }
}

void game_unload(struct game *game)
{
    if (game->handle) {
        game->api.finalize(game->state);
        game->state = NULL;
        dlclose(game->handle);
        game->handle = NULL;
        game->id = 0;
    }
}

void* hot_loop_batch(std::string arg_so_location) {
  strcpy(SO_LOCATION, arg_so_location.c_str());
  struct game game = {0};
  if (SO_LOCATION[0]) {
    game_load(&game, false /* bInteractive */);
  }
  else {
    printf("In hot_loop_batch: handle is not valid --> exit;\n");
    exit(1);
  }
  if (game.handle)
    game.api.step(game.state);
  else {
    printf("In hot_loop_batch: handle is not valid --> exit;\n");
    exit(1);
  }
  return game.state;
}

void* hot_loop_interactive(std::string arg_so_location) {
  strcpy(SO_LOCATION, arg_so_location.c_str());
  struct game game = {0};

  long dt;
  long started;
  int running = 1;
  XWindow xw;
  struct nk_context *ctx;

  /* X11 */
  memset(&xw, 0, sizeof xw);
  xw.dpy = XOpenDisplay(NULL);
  if (!xw.dpy) die("Could not open a display; perhaps $DISPLAY is not set?");
  xw.root = DefaultRootWindow(xw.dpy);
  xw.screen = XDefaultScreen(xw.dpy);
  xw.vis = XDefaultVisual(xw.dpy, xw.screen);
  xw.cmap = XCreateColormap(xw.dpy,xw.root,xw.vis,AllocNone);

  xw.swa.colormap = xw.cmap;
  xw.swa.event_mask =
    ExposureMask | KeyPressMask | KeyReleaseMask |
    ButtonPress | ButtonReleaseMask| ButtonMotionMask |
    Button1MotionMask | Button3MotionMask | Button4MotionMask | Button5MotionMask|
    PointerMotionMask | KeymapStateMask;
  xw.win = XCreateWindow(xw.dpy, xw.root, 0, 0, WINDOW_WIDTH, WINDOW_HEIGHT, 0,
                         XDefaultDepth(xw.dpy, xw.screen), InputOutput,
                         xw.vis, CWEventMask | CWColormap, &xw.swa);

  XStoreName(xw.dpy, xw.win, "X11");
  XMapWindow(xw.dpy, xw.win);
  xw.wm_delete_window = XInternAtom(xw.dpy, "WM_DELETE_WINDOW", False);
  XSetWMProtocols(xw.dpy, xw.win, &xw.wm_delete_window, 1);
  XGetWindowAttributes(xw.dpy, xw.win, &xw.attr);
  xw.width = (unsigned int)xw.attr.width;
  xw.height = (unsigned int)xw.attr.height;

  /* GUI */
  xw.font = nk_xfont_create(xw.dpy, "fixed");
  ctx = nk_xlib_init(xw.font, xw.dpy, xw.screen, xw.win, xw.width, xw.height);

  /* style.c */
  /*set_style(ctx, THEME_WHITE);*/
  /* set_style(ctx, THEME_RED); */
  /*set_style(ctx, THEME_BLUE);*/
  /*set_style(ctx, THEME_DARK);*/

  while (running) {
    /* Input */
    XEvent evt;
    started = timestamp();
    nk_input_begin(ctx);
    while (XPending(xw.dpy)) {
      XNextEvent(xw.dpy, &evt);
      if (evt.type == ClientMessage) goto cleanup;
      if (XFilterEvent(&evt, xw.win)) continue;
      nk_xlib_handle_event(xw.dpy, xw.screen, xw.win, &evt);
    }
    nk_input_end(ctx);

    /* GUI */
    // Empty: moved to game_step.
    if (SO_LOCATION[0]) {
      printf("I would like to load %s\n", SO_LOCATION);
      game_load(&game, true /* bInteractive */);
      *SO_LOCATION = 0;
    }
    if (game.handle)
      if (!game.api.step(game.state))
        break;

    if (nk_window_is_hidden(ctx, "Plot Controls")) break;

    /* Draw */
    XClearWindow(xw.dpy, xw.win);
    nk_xlib_render(xw.win, nk_rgb(30,30,30));
    XFlush(xw.dpy);


    /* Timing */
    dt = timestamp() - started;
    if (dt < DTIME)
      sleep_for(DTIME - dt);
  }

 cleanup:
  nk_xfont_del(xw.dpy, xw.font);
  nk_xlib_shutdown();
  XUnmapWindow(xw.dpy, xw.win);
  XFreeColormap(xw.dpy, xw.cmap);
  XDestroyWindow(xw.dpy, xw.win);
  XCloseDisplay(xw.dpy);
  return game.state;
}

void *hot_loop(std::string arg_so_location, bool bInteractive) {
  void *result = NULL;
  if (bInteractive)
    result = hot_loop_interactive(arg_so_location);
  else
    result = hot_loop_batch(arg_so_location);
  printf("In hot_loop...\n");
  return result;
}
