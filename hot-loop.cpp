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

char SO_LOCATION[256] = "";

char *errstr;

struct game {
    void *handle;
    ino_t id;
    struct game_api api;
    struct game_state *state;
};

static void game_load(struct game *game)
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
        // struct game_api *api = &GAME_API;
        if (api != NULL) {
          game->api = *api;
          if (game->state == NULL) game->state = game->api.init();
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

int hot_loop(std::string arg_so_location)
{
  strcpy(SO_LOCATION, arg_so_location.c_str());
  // TApplication *theApp = new TApplication("tapp", NULL, NULL);
  // theApp = theApp;
  struct game game = {0};
  for (;;) {
    if (SO_LOCATION[0]) {
      printf("I would like to load %s\n", SO_LOCATION);
      game_load(&game);
      *SO_LOCATION = 0;
    }
    if (game.handle)
      if (!game.api.step(game.state))
        break;
    usleep(200000);
  }
  game_unload(&game);
  return 0;
}
