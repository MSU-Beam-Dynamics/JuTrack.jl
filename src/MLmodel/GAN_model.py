"""
GAN model for space-charge Hamioltonian generation.
Related paper: "A symplectic machine learning model for fast simulation of space-charge effects" by J. Wan, Y. Hao and J. Qiang.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np

# Define the Generator using U-Net architecture
class GeneratorUNet(nn.Module):
    def __init__(self, in_channels=1, out_channels=1):
        super(GeneratorUNet, self).__init__()

        # Encoder (Downsampling layers)
        self.down1 = self.down_block(in_channels, 64, normalize=False)  # 128 -> 64
        self.down2 = self.down_block(64, 128)                           # 64 -> 32
        self.down3 = self.down_block(128, 256)                          # 32 -> 16
        self.down4 = self.down_block(256, 512)                          # 16 -> 8
        self.down5 = self.down_block(512, 512)                          # 8 -> 4
        self.down6 = self.down_block(512, 512)                          # 4 -> 2
        self.down7 = self.down_block(512, 512, normalize=False)         # 2 -> 1

        # Decoder (Upsampling layers)
        self.up1 = self.up_block(512, 512, dropout=0.5)                 # 1 -> 2
        self.up2 = self.up_block(1024, 512, dropout=0.5)                # Input channels after concatenation
        self.up3 = self.up_block(1024, 512, dropout=0.5)
        self.up4 = self.up_block(1024, 512, dropout=0.5)
        self.up5 = self.up_block(768, 256)                              
        self.up6 = self.up_block(384, 128)                              
        self.up7 = self.up_block(192, 64)                               

        # Final output layer
        self.final = nn.Sequential(
            nn.Conv2d(64, out_channels, kernel_size=3, padding=1)
        )

    def down_block(self, in_channels, out_channels, normalize=True):
        """Defines a downsampling block"""
        layers = [nn.Conv2d(in_channels, out_channels, kernel_size=4, stride=2, padding=1)]
        if normalize:
            layers.append(nn.InstanceNorm2d(out_channels))
        layers.append(nn.LeakyReLU(0.2, inplace=True))
        return nn.Sequential(*layers)

    def up_block(self, in_channels, out_channels, dropout=0.0):
        """Defines an upsampling block"""
        layers = [
            nn.ConvTranspose2d(in_channels, out_channels, kernel_size=4, stride=2, padding=1),
            nn.InstanceNorm2d(out_channels),
            nn.ReLU(inplace=True)
        ]
        if dropout:
            layers.append(nn.Dropout(dropout))
        return nn.Sequential(*layers)

    def forward(self, x):
        # Encoder
        d1 = self.down1(x)   # 128 -> 64
        d2 = self.down2(d1)  # 64 -> 32
        d3 = self.down3(d2)  # 32 -> 16
        d4 = self.down4(d3)  # 16 -> 8
        d5 = self.down5(d4)  # 8 -> 4
        d6 = self.down6(d5)  # 4 -> 2
        d7 = self.down7(d6)  # 2 -> 1

        # Decoder with skip connections
        u1 = self.up1(d7)                  # 1 -> 2
        u1 = torch.cat([u1, d6], 1)        # Channels: 512 + 512 = 1024

        u2 = self.up2(u1)                  # 2 -> 4
        u2 = torch.cat([u2, d5], 1)        # Channels: 512 + 512 = 1024

        u3 = self.up3(u2)                  # 4 -> 8
        u3 = torch.cat([u3, d4], 1)        # Channels: 512 + 512 = 1024

        u4 = self.up4(u3)                  # 8 -> 16
        u4 = torch.cat([u4, d3], 1)        # Channels: 512 + 256 = 768

        u5 = self.up5(u4)                  # 16 -> 32
        u5 = torch.cat([u5, d2], 1)        # Channels: 256 + 128 = 384

        u6 = self.up6(u5)                  # 32 -> 64
        u6 = torch.cat([u6, d1], 1)        # Channels: 128 + 64 = 192

        u7 = self.up7(u6)                  # 64 -> 128

        # Final output
        output = self.final(u7)            # Output size: [batch_size, out_channels, 128, 128]
        return output


# Define the Discriminator using a PatchGAN
class Discriminator(nn.Module):
    def __init__(self, in_channels=1):
        super(Discriminator, self).__init__()

        def discriminator_block(in_filters, out_filters, normalize=True):
            """Defines a discriminator block"""
            layers = [nn.Conv2d(in_filters, out_filters, kernel_size=4, stride=2, padding=1)]
            if normalize:
                layers.append(nn.InstanceNorm2d(out_filters))
            layers.append(nn.LeakyReLU(0.2, inplace=True))
            return layers

        # PatchGAN discriminator layers
        self.model = nn.Sequential(
            *discriminator_block(in_channels * 2, 64, normalize=False),
            *discriminator_block(64, 128),
            *discriminator_block(128, 256),
            # *discriminator_block(256, 512),
            nn.ZeroPad2d((1, 0, 1, 0)),  
            nn.Conv2d(256, 1, kernel_size=4, padding=1)  
        )

    def forward(self, img_A, img_B):
        # Concatenate input and output images by channels
        img_input = torch.cat((img_A, img_B), 1)
        return self.model(img_input)

def smooth_prediction(prediction, kernel_size=5, sigma=2):
    # Smooth the prediction using a Gaussian kernel
    prediction = torch.FloatTensor(prediction)
    coords = torch.arange(kernel_size, dtype=torch.float32) - (kernel_size - 1) / 2
    grid = torch.exp(-0.5 * (coords ** 2) / (sigma ** 2))
    kernel_1d = grid / grid.sum()
    kernel_2d = kernel_1d[:, None] @ kernel_1d[None, :]
    kernel_2d = kernel_2d.to(prediction.device).view(1, 1, kernel_size, kernel_size)
    
    smoothed_prediction = F.conv2d(prediction, kernel_2d, padding=kernel_size // 2)
    return smoothed_prediction.numpy()
