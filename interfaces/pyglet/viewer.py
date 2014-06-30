import pyglet

from src import map_generator
import map_view


class YasgWindow(pyglet.window.Window):
    def __init__(self):
        super(YasgWindow, self).__init__(800, 600, "Yasg")
        
        self.world = map_generator.generate_map(1000)
        self.map_view = map_view.MapView(self.world)
        
        self.push_handlers(self.map_view)

    def on_draw(self):
        self.clear()
        
        self.map_view.draw()

if __name__ == '__main__':
    game = YasgWindow()
    pyglet.app.run()