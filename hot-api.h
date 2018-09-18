#pragma once

#include <stdbool.h>

struct game_state;

struct game_api {
    /**
     * @return a fresh game state
     */
    void *(*init)();

    /**
     * Destroys a game state.
     */
    void (*finalize)(void *state);

    /**
     * Called exactly once when the game code is reloaded.
     */
    void (*reload)(void *state);

    /**
     * Called exactly once when the game code is about to be reloaded.
     */
    void (*unload)(void *state);

    /**
     * Called at a regular interval by the main program.
     * @return true if the program should continue
     */
    bool (*step)(void *state);
};

extern const struct game_api GAME_API;
// extern struct game_api GAME_API;
