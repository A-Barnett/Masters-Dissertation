#include "ProceduralTerrain.h"
#include "DrawDebugHelpers.h"
#include "Misc/DateTime.h"
#include "Misc/Timespan.h"
#include "Containers/UnrealString.h"
#include "ProceduralMeshComponent.h"
#include "Async/Async.h"
#include "Engine/World.h"            
#include "PhysicsEngine/BodySetup.h" 
#include "PhysicsEngine/ConvexElem.h" 
#include "Containers/Array.h"       
#include "Math/Vector.h"           
#include "Math/Vector2D.h"           
#include "Math/Color.h"            
#include "Components/SceneComponent.h"
#include "GameFramework/Character.h"
#include "GameFramework/Pawn.h"
#include "GameFramework/PlayerController.h"
#include "GameFramework/Actor.h"
#include "GameFramework/Controller.h"
#include "GameFramework/MovementComponent.h"
#include "Components/InputComponent.h"
#include "GameFramework/CharacterMovementComponent.h"


AProceduralTerrain::AProceduralTerrain()
{
	PrimaryActorTick.bCanEverTick = true;
}

// Called when the game starts or when spawned
void AProceduralTerrain::BeginPlay()
{
	Super::BeginPlay();
	randomOffset.Value = rand() % 100000;
	randomOffset.Key = rand() % 100000;
	MeshComponent.Empty();
	TerrainComponents.Empty();
	PlayerPawn = GetWorld()->GetFirstPlayerController()->GetPawn();
	currentGridSize = (Width+1) * Scale;
	PlayerGridPos = FVector2D(0, 0);
	RealtimeMesh = GetRealtimeMeshComponent()->InitializeRealtimeMesh<URealtimeMeshSimple>();
	RealtimeMesh->SetupMaterialSlot(0, "GrassMat", TerrainMaterial);
	RealtimeMesh->SetupMaterialSlot(1, "PathMat", PathMaterial);
	FRealtimeMeshCollisionConfiguration CollisionConfig;
	CollisionConfig.bUseAsyncCook = true;  // Allows async cooking for better performance
	CollisionConfig.bUseComplexAsSimpleCollision = true;
	CollisionConfig.bShouldFastCookMeshes = true;
	RealtimeMeshComponent->SetMobility(EComponentMobility::Static);
	RealtimeMeshComponent->SetCastShadow(false);
	RealtimeMeshComponent->SetGenerateOverlapEvents(false);
	RealtimeMeshComponent->SetCollisionProfileName(UCollisionProfile::BlockAll_ProfileName);
	RealtimeMeshComponent->SetCollisionEnabled(ECollisionEnabled::QueryAndPhysics);
	SetRootComponent(RealtimeMeshComponent);
	RealtimeMesh->SetCollisionConfig(CollisionConfig);


	GenerateTerrain();

}



void AProceduralTerrain::GenerateTerrain()
{
	FDateTime StartTime = FDateTime::Now();
	GeneratePath(true,0);
	SmoothPathPointsHeight(PathHeightSmooth,0);
	GeneratePathMesh(0);
	pointsGenerated += NumPoints;
	UpdateTerrain(true);

	FDateTime EndTime = FDateTime::Now();
	FTimespan Duration = EndTime - StartTime;
	UE_LOG(LogTemp, Display, TEXT("TERRAIN GEN COMPLETE IN %f SECONDS"), Duration.GetTotalSeconds());
}

// Called every frame
void AProceduralTerrain::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);
	if (GetWorld()->GetFirstPlayerController()->IsInputKeyDown(EKeys::M)) {
		ResetPlayerPhysics();
	}

	FVector PlayerPos = PlayerPawn->GetActorLocation();
	if (PlayerPos.X > ((PlayerGridPos.X * currentGridSize) + (currentGridSize))) {
		//right
		PlayerGridPos.X++;
		UpdateTerrain();
		UE_LOG(LogTemp, Display, TEXT("RIGHT"));
	}
	else if (PlayerPos.X < ((PlayerGridPos.X * currentGridSize))) {
		//left
		PlayerGridPos.X--;
		UpdateTerrain();
		UE_LOG(LogTemp, Display, TEXT("LEFT"));
	}
	else if (PlayerPos.Y > ((PlayerGridPos.Y * currentGridSize) + currentGridSize)) {
		//top
		PlayerGridPos.Y++;
		UpdateTerrain();
		UE_LOG(LogTemp, Display, TEXT("TOP"));
	}
	else if (PlayerPos.Y < ((PlayerGridPos.Y * currentGridSize))) {
		//bottom
		PlayerGridPos.Y--;
		UpdateTerrain();
		UE_LOG(LogTemp, Display, TEXT("BOTTOM"));
	}

}

void AProceduralTerrain::ResetPlayerPhysics()
{
	if (!PlayerPawn) return;

	// Reset rotation with Teleport flag
	PlayerPawn->SetActorRotation(FRotator::ZeroRotator, ETeleportType::TeleportPhysics);

	if (ACharacter* Character = Cast<ACharacter>(PlayerPawn))
	{
		if (UCharacterMovementComponent* MovementComp = Character->GetCharacterMovement())
		{
			MovementComp->Velocity = FVector::ZeroVector;
		}
	}

	if (UPrimitiveComponent* RootPrimitive = Cast<UPrimitiveComponent>(PlayerPawn->GetRootComponent()))
	{
		RootPrimitive->SetPhysicsLinearVelocity(FVector::ZeroVector);
		RootPrimitive->SetPhysicsAngularVelocityInDegrees(FVector::ZeroVector);
	}
}




void AProceduralTerrain::UpdateTerrain(bool initial)
{
	if (!initial) {
		GeneratePath(false, pointsGenerated);
		SmoothPathPointsHeight(PathHeightSmooth, pointsGenerated);
		GeneratePathMesh(pointsGenerated);

		pointsGenerated += makePerUpdate;
	}
	UE_LOG(LogTemp, Display, TEXT("POINTS %i"), PathPoints.Num());
	TArray<TerrainComponent*> terrainsToBeMoved;
	TArray<TerrainComponent*> ComponentsToUpdate;
	TArray<std::pair<FVector2D, int>> terrainsToMake;

	const float range4 = 16.0f;
	const int distanceSearch = range4 + 4;

	// Define LOD ranges: {distanceLimit, LOD}
	const TArray<std::pair<float, int>> LODRanges = {
		{2.0f, 0},
		{4.0f, 1},
		{7.0f, 2},
		{10.0f, 3},
		{range4, 4}
	};

	auto FindLODForDistance = [&](float distance) -> int
		{
			for (const auto& pair : LODRanges)
			{
				if (distance <= pair.first)
					return pair.second;
			}
			return -1; // Out of range
		};

	for (int x = -distanceSearch; x <= distanceSearch; ++x)
	{
		for (int y = -distanceSearch; y <= distanceSearch; ++y)
		{
			FVector2D GridPosition = PlayerGridPos + FVector2D(x, y);
			const float distance = FVector2D::Distance(GridPosition, PlayerGridPos);
			const int targetLOD = FindLODForDistance(distance);

			TerrainComponent* component = FindTerrainComponent(GridPosition);

			if (targetLOD != -1)
			{
				if (component)
				{
					if (component->GetLOD() != targetLOD)
					{
						component->SetLOD(targetLOD);
						component->SetGridPosition(GridPosition);
						component->SetIsNewPos(false);
						ComponentsToUpdate.Add(component);
					}
				}
				else
				{
					terrainsToMake.Add({ GridPosition, targetLOD });
				}
			}
			else if (component)
			{
				terrainsToBeMoved.Add(component);
			}
		}
	}

	// Try to reuse old terrain components before creating new ones
	int terrainPiecesMade = 0;
	for (int i = 0; i < terrainsToMake.Num(); ++i)
	{
		const FVector2D& GridPos = terrainsToMake[i].first;
		const int LOD = terrainsToMake[i].second;

		if (i < terrainsToBeMoved.Num())
		{
			TerrainComponent* reusable = terrainsToBeMoved[i];
			reusable->SetGridPosition(GridPos);
			reusable->SetLOD(LOD);
			reusable->SetIsNewPos(true);
			ComponentsToUpdate.Add(reusable);
		}
		else
		{
			ComponentsToUpdate.Add(CreateTerrainComponent(GridPos, LOD));
			++terrainPiecesMade;
		}
	}

	// Logging
	UE_LOG(LogTemp, Display, TEXT("New Terrains needed: %i"), terrainsToMake.Num());
	UE_LOG(LogTemp, Display, TEXT("Old Terrains no longer needed: %i"), terrainsToBeMoved.Num());
	UE_LOG(LogTemp, Display, TEXT("New Terrains made: %i"), terrainPiecesMade);
	UE_LOG(LogTemp, Display, TEXT("Total Terrains updated: %i"), ComponentsToUpdate.Num());
	ParallelFor(ComponentsToUpdate.Num(), [this, ComponentsToUpdate](int32 Index)
	{
		GenerateTerrainSection(ComponentsToUpdate[Index]);
	});
}


TerrainComponent* AProceduralTerrain::FindTerrainComponent(const FVector2D& GridPosition) {
	for (TerrainComponent* Component : TerrainComponents)
	{
		if (Component && Component->GetGridPosition() == GridPosition)
		{
			return Component;
		}
	}
	return nullptr;
}

TerrainComponent* AProceduralTerrain::CreateTerrainComponent(const FVector2D& GridPosition, int LOD)
{
	totalCreated++;
	FString ComponentName = FString::Printf(TEXT("MeshComponent%d"), TerrainComponents.Num());
	TerrainComponent* NewComponent = new TerrainComponent(nullptr, GridPosition, LOD, TerrainComponents.Num()); //nullptr
	NewComponent->SetIsActive(true);
	TerrainComponents.Add(NewComponent);
	return NewComponent;
}


float AProceduralTerrain::CalculateNoiseAtPoint(int32 X, int32 Y) const {
	float FinalHeight = 0.0f;
	for (const FTerrainLayer& Layer : TerrainLayers)
	{
		if (Layer.Enabled) {
			float LayerHeight = FMath::PerlinNoise2D(FVector2D(X + randomOffset.Value, Y + randomOffset.Key) * (Layer.Frequency / 4)) * Layer.Amplitude;
			FinalHeight += LayerHeight * Layer.Contribution;
		}
	}
	return FinalHeight;
}


float AProceduralTerrain::CalculateHeightOnPath(int32 X, int32 Y, TArray<FVector3f>* VertMap) const
{

	float closestDist = INFINITY;
	FVector3f closestPoint;
	FVector3f p = FVector3f(X * Scale, Y * Scale, 0);
	bool PointFound = false;
	for (FVector3f point : *VertMap) {
		if (FVector3f::DistXY(p, point) < closestDist) {
			closestDist = FVector3f::Dist2D(p, point);
			closestPoint = point;
		}
	}
	return closestPoint.Z - (HeightAdjust * 2);
}


float AProceduralTerrain::CalculateHeight(int32 X, int32 Y, TArray<FVector3f>* PointsMap, TArray<FVector3f>* VertMap) const
{
	float pointHeight = CalculateNoiseAtPoint(X, Y);
	if (PointsMap == nullptr) return pointHeight;
	if (IsOnPath(X, Y, true, PointsMap))
	{
		return CalculateHeightOnPath(X, Y, VertMap);
	}
	else if (IsOnPath(X, Y, false, PointsMap)) {
		//return AverageHeightNeighbours(X, Y, SmoothingSize);
		float pathHeight = CalculateHeightOnPath(X, Y, VertMap);
		float lerpAmount = (DistFromPath(X, Y, false, PointsMap) / (SmoothingThicknessOffset - ThicknessOffset));
		return FMath::Lerp(pathHeight, pointHeight, lerpAmount);
	}
	return pointHeight;
}


bool AProceduralTerrain::IsOnPath(int32 X, int32 Y, bool useOffset, TArray<FVector3f>* PointsMap) const
{
	FVector2D Point2D(X * Scale, Y * Scale);
	FVector Point(Point2D.X, Point2D.Y, 0.0f);
	for (int32 i = 0; i < PointsMap->Num()-1; ++i)
	{
		// Get current and next point, converting to FVector
		FVector2D CurrentPoint2D = FVector2D((*PointsMap)[i].X, (*PointsMap)[i].Y);
		FVector2D NextPoint2D = FVector2D((*PointsMap)[(i + 1) ].X, (*PointsMap)[(i + 1) ].Y);

		FVector CurrentPoint(CurrentPoint2D.X, CurrentPoint2D.Y, 0.0f);
		FVector NextPoint(NextPoint2D.X, NextPoint2D.Y, 0.0f);

		// Calculate the distance from the point to the segment
		float Distance = FMath::PointDistToSegment(Point, CurrentPoint, NextPoint);

		if (Distance <= (useOffset ? (Thickness + ThicknessOffset) : (Thickness + SmoothingThicknessOffset)))
		{
			return true;
		}
	}
	return false; // Point is not near the path
}


float AProceduralTerrain::DistFromPath(int32 X, int32 Y, bool useOffset, TArray<FVector3f>* PointsMap) const
{
	FVector2D Point2D(X * Scale, Y * Scale);
	FVector Point(Point2D.X, Point2D.Y, 0.0f);
	float FinalDistance = INFINITY;
	for (int32 i = 0; i < PointsMap->Num()-1; ++i)
	{
		// Get current and next point, converting to FVector
		FVector2D CurrentPoint2D = FVector2D((*PointsMap)[i].X, (*PointsMap)[i].Y);
		FVector2D NextPoint2D = FVector2D((*PointsMap)[(i + 1)].X, (*PointsMap)[(i + 1)].Y);

		FVector CurrentPoint(CurrentPoint2D.X, CurrentPoint2D.Y, 0.0f);
		FVector NextPoint(NextPoint2D.X, NextPoint2D.Y, 0.0f);

		// Calculate the distance from the point to the segment
		float Distance = FMath::PointDistToSegment(Point, CurrentPoint, NextPoint);

		if (Distance <= FinalDistance)
		{
			FinalDistance = Distance;
		}
	}

	return (FinalDistance - (Thickness + ThicknessOffset));

}


FVector2D CatmullRomInterpolate(const FVector2D& P0, const FVector2D& P1, const FVector2D& P2, const FVector2D& P3, float T)
{
	float T2 = T * T;
	float T3 = T2 * T;

	return 0.5f * (
		(2.0f * P1) +
		(-P0 + P2) * T +
		(2.0f * P0 - 5.0f * P1 + 4.0f * P2 - P3) * T2 +
		(-P0 + 3.0f * P1 - 3.0f * P2 + P3) * T3
		);
}


void AProceduralTerrain::GeneratePath(bool start, int32 offset)
{
	// PathPoints.Empty();

	TArray<FVector2D> BasePoints;
	float StepSize = Scale * 5.0f; // Step distance per point
	float MaxTurnAngle = 2.5f;
	int32 pointsToMake = makePerUpdate+1;

	if (start) {
		RandomStream = FRandomStream(PathSeed);
		CurrentPosition = FVector2D(-Width * Scale * 20.0f, 0.0f); // Start position
		CurrentDirection = FVector2D(1.0f, 0.0f); // Initial direction (X-axis)
		CurrentTurnAngle = 0.0f;
		BasePoints.Add(CurrentPosition);
		pointsToMake = NumPoints;
	}

	for (int32 i = offset+1; i < pointsToMake+ offset; ++i)
	{
		if (i % TurnSize == 0)
		{
			float NewTurnAngle = RandomStream.FRandRange(-MaxTurnAngle, MaxTurnAngle);
			CurrentTurnAngle = NewTurnAngle / TurnSmoothPoints; // Spread it over multiple points
		}

		// Apply the gradual turn
		float AngleRad = FMath::DegreesToRadians(CurrentTurnAngle);
		float CosA = FMath::Cos(AngleRad);
		float SinA = FMath::Sin(AngleRad);

		FVector2D NewDirection = FVector2D(
			CurrentDirection.X * CosA - CurrentDirection.Y * SinA,
			CurrentDirection.X * SinA + CurrentDirection.Y * CosA
		).GetSafeNormal();

		CurrentDirection = NewDirection;
		CurrentPosition += CurrentDirection * StepSize;
		BasePoints.Add(CurrentPosition);
	}

	// Apply Catmull-Rom interpolation to smooth the path further
	for (int32 i = 0; i < BasePoints.Num(); ++i)
	{
		const FVector2D& P0 = BasePoints[FMath::Max(0, i - 1)];
		const FVector2D& P1 = BasePoints[i];
		const FVector2D& P2 = BasePoints[FMath::Min(i + 1, BasePoints.Num() - 1)];
		const FVector2D& P3 = BasePoints[FMath::Min(i + 2, BasePoints.Num() - 1)];

		float T = 1.0f;
		FVector2D SmoothedPoint = CatmullRomInterpolate(P0, P1, P2, P3, T);
		FVector3f FinalPoint = FVector3f(SmoothedPoint.X, SmoothedPoint.Y, CalculateNoiseAtPoint(SmoothedPoint.X / Scale, SmoothedPoint.Y / Scale));
		
		PathPoints.Add(FinalPoint);
		int32 GridX = FMath::FloorToInt(FinalPoint.X / currentGridSize);
		int32 GridY = FMath::FloorToInt(FinalPoint.Y / currentGridSize);
		UE_LOG(LogTemp, Display, TEXT("Point at Grid (%i,%i)"), GridX, GridY);
		FIntPoint Cell = FIntPoint(GridX, GridY);
		PathGridMap.FindOrAdd(Cell).Add(FinalPoint);
			
	}
}


void AProceduralTerrain::SmoothPathPointsHeight(float smoothLevel, int32 offset)
{
	if (offset != 0) offset--;
	UE_LOG(LogTemp, Display, TEXT("Offset %i, Points %i"),offset,PathPoints.Num());
	for (int32 i = 1+ offset; i < (PathPoints.Num()); ++i)
	{
		FVector3f& CurrentPoint = PathPoints[i];
		FVector3f& PrevPoint = PathPoints[i - 1];
		float InterpolatedHeight = FMath::Lerp(CurrentPoint.Z, PrevPoint.Z, smoothLevel);
		CurrentPoint.Z = InterpolatedHeight;
	}
}


void AProceduralTerrain::GenerateTerrainSection(TerrainComponent* Component)
{
	// Launch async task for heavy computation
	Async(EAsyncExecution::ThreadPool, [=, this]()
		{

			////////////////////////////////////
			TSharedPtr<TRealtimeMeshBuilderLocal<uint16, FPackedNormal, FVector2DHalf, 1>> StoredBuilder = Component->GetPrevBuilder();

			TSharedPtr<FRealtimeMeshStreamSet> StreamSet = MakeShared<FRealtimeMeshStreamSet>();
			auto BuilderPtr = MakeShared<TRealtimeMeshBuilderLocal<uint16, FPackedNormal, FVector2DHalf, 1>>(*StreamSet);
			auto& Builder = *BuilderPtr;;
		
			// here we go ahead and enable all the basic mesh data parts
			Builder.EnableTangents();
			Builder.EnableTexCoords();
			Builder.EnablePolyGroups();


			int32 SectionSize = Width;
			int32 SectionX = Component->GetGridPosition().X;
			int32 SectionY = Component->GetGridPosition().Y;
			int32 StartX = SectionX * SectionSize;
			int32 StartY = SectionY * SectionSize;
			int32 EndX = StartX + SectionSize;
			int32 EndY = StartY + SectionSize;
			int32 LOD = pow(2, Component->GetLOD());
			int32 vertsPerRow = (SectionSize / LOD) + 1;
			int32 totalSize = vertsPerRow * vertsPerRow;

			FIntPoint QueryCell = FIntPoint(Component->GetGridPosition().X,Component->GetGridPosition().Y);
			TArray<FVector3f>* PointsMap = new TArray<FVector3f>();
			TArray<FVector3f>* VertsMap = new TArray<FVector3f>();
			for (int32 dx = -2; dx <= 2; ++dx)
			{
				for (int32 dy = -2; dy <= 2; ++dy)
				{
					FIntPoint NeighborCell = FIntPoint(QueryCell.X + dx, QueryCell.Y + dy);
					if (TArray<FVector3f>* NeighborPoints = PathGridMap.Find(NeighborCell))
					{
						PointsMap->Append(*NeighborPoints);
					}
					if (TArray<FVector3f>* NeighborVerts = PathVertMap.Find(NeighborCell))
					{
						VertsMap->Append(*NeighborVerts);
					}
				}
			}

			// Generate vertices and UVs
			int vertsAmount = 0;
			if (Component->TryLockMesh() && StoredBuilder.Get() && StoredBuilder.Get()->NumVertices() > 0 && !Component->GetIsNewPos()) {
				int32 OldLOD = pow(2, Component->GetOldLOD());
				int32 pointsOutsideRange = 0;
				if (Component->GetOldLOD() < Component->GetLOD()) {
					int32 SamplingFactor = LOD / OldLOD; // Downsample ratio
					int32 OldvertsPerRow = (SectionSize / OldLOD) + 1;
					for (int32 y = 0; y < OldvertsPerRow; y += 2) {
						for (int32 x = 0; x < OldvertsPerRow; x += 2) {
							int32 OldIndex = y * OldvertsPerRow + x; // Map old vertex to reduced grid
							if (StoredBuilder && vertsAmount < StoredBuilder.Get()->NumVertices()) {
								Builder.AddVertex(StoredBuilder.Get()->GetPosition(OldIndex)).SetTexCoord(FVector2f(
									static_cast<float>(x) / SectionSize,
									static_cast<float>(y) / SectionSize
								) * UVScale);
							}
							else {
								float Z = CalculateHeight(x, y, PointsMap, VertsMap);
								Builder.AddVertex(FVector3f(x * Scale, y * Scale, Z))
									.SetTexCoord(FVector2f(
										static_cast<float>(x/ SectionSize),
										static_cast<float>(y  / SectionSize)
									) * UVScale);
								pointsOutsideRange++;
							}
							vertsAmount++;
						}
					}
					if (pointsOutsideRange > 0) {
						UE_LOG(LogTemp, Display, TEXT("OLD LOD %i, NEW LOD %i ,Points Failed %i ,Total Points %i, Sample Factor %i"), Component->GetOldLOD(), Component->GetLOD(), pointsOutsideRange, vertsAmount, SamplingFactor);
					}
				}
				else if(Component->GetOldLOD() > Component->GetLOD()) {
					int32 SamplingFactor = OldLOD / LOD; // Upsample ratio
					int32 OldvertsPerRow = (SectionSize / OldLOD) + 1;
					int32 OldtotalSize = OldvertsPerRow * OldvertsPerRow; // Total size of stored array
					int32 pointsUsed = 0;
			
					for (int32 y = StartY; y <= EndY; y += LOD) {
						for (int32 x = StartX; x <= EndX; x += LOD) {
							// Map to old vertex position if available
							if (StoredBuilder && y % OldLOD == 0 && x % OldLOD == 0 && pointsUsed < StoredBuilder.Get()->NumVertices()) {
								// Safe index check and sampling
								Builder.AddVertex(StoredBuilder.Get()->GetPosition(pointsUsed)).SetTexCoord(FVector2f(
									static_cast<float>(x - StartX) / SectionSize,
									static_cast<float>(y - StartY) / SectionSize
								) * UVScale); // Sample from old vertex
								pointsUsed++;
							} else {
								if (x % OldLOD == 0 && y % OldLOD == 0) {
									pointsOutsideRange++;
								}
								// Generate new vertex if out of bounds or misaligned
								float Z = CalculateHeight(x, y, PointsMap, VertsMap);
								Builder.AddVertex(FVector3f(x * Scale, y * Scale, Z))
									.SetTexCoord(FVector2f(
										static_cast<float>(x - StartX) / SectionSize,
										static_cast<float>(y - StartY) / SectionSize
									) * UVScale);
							}
							vertsAmount++;
						}
					}
					if (pointsOutsideRange > 0) {
						UE_LOG(LogTemp, Display, TEXT("OLD LOD %i, NEW LOD %i ,Points Used %i ,Points Failed %i ,Total Points %i, Sample Factor %i"), Component->GetOldLOD(), Component->GetLOD(), pointsUsed, pointsOutsideRange, vertsAmount, OldLOD);
					}
					
				}
				Component->UnlockMesh();
				
			}
			else {
				for (int32 y = StartY; y <= EndY; y += LOD)
				{
					for (int32 x = StartX; x <= EndX; x += LOD)
					{
						float Z = CalculateHeight(x, y, PointsMap, VertsMap);
						Builder.AddVertex(FVector3f(x * Scale, y * Scale, Z))
							.SetTexCoord(FVector2f(
								static_cast<float>(x - StartX) / SectionSize,
								static_cast<float>(y - StartY) / SectionSize
							) * UVScale);
						vertsAmount++;
					}
				}
			}

			int triangleCount = 0;
			for (int32 y = 0; y < SectionSize; y += LOD)
			{
				for (int32 x = 0; x < SectionSize; x += LOD)
				{
					int32 CurrentRow = y / LOD;
					int32 NextRow = (y + LOD) / LOD;

					int32 CurrentIndex = CurrentRow * vertsPerRow + (x / LOD);
					int32 RightIndex = CurrentIndex + 1;
					int32 BottomIndex = NextRow * vertsPerRow + (x / LOD);
					int32 BottomRightIndex = BottomIndex + 1;

					Builder.AddTriangle(CurrentIndex, BottomIndex, RightIndex, 0);
					Builder.AddTriangle(RightIndex, BottomIndex, BottomRightIndex, 0);
					triangleCount += 2;
				}
			}


			// Calculate normals and tangents without using texture coordinates
			TArray<FVector3f> VertexNormals;
			VertexNormals.SetNum(totalSize, false);
			TArray<FVector3f> VertexTangents;
			VertexTangents.SetNum(totalSize, false);
			for (int32 i = 0; i < totalSize; i++) {
				VertexNormals[i] = FVector3f::ZeroVector;
				VertexTangents[i] = FVector3f::ZeroVector;
			}
			if (Builder.NumVertices() <= 0) {
				return;
			}

			for (int32 i = 0; i < triangleCount; i++)
			{

				TIndex3<uint32> Index = Builder.GetTriangle(i);
				const FVector3f& Vertex0 = Builder.GetPosition(Index[0]);
				const FVector3f& Vertex1 = Builder.GetPosition(Index[1]);
				const FVector3f& Vertex2 = Builder.GetPosition(Index[2]);

				// Calculate edges
				FVector3f Edge1 = Vertex1 - Vertex0;
				FVector3f Edge2 = Vertex2 - Vertex0;

				// Calculate the face normal using the cross product
				FVector3f FaceNormal = FVector3f::CrossProduct(Edge1, Edge2).GetSafeNormal() * -1.0f;

				// Accumulate normals for each vertex in the triangle
				VertexNormals[Index[0]] += FaceNormal;
				VertexNormals[Index[1]] += FaceNormal;
				VertexNormals[Index[2]] += FaceNormal;

				// Generate a tangent using an arbitrary orthogonal vector to the normal
				FVector3f Tangent = FVector3f::CrossProduct(FaceNormal, Edge1).GetSafeNormal();

				// Accumulate tangents for each vertex in the triangle
				VertexTangents[Index[0]] += Tangent;
				VertexTangents[Index[1]] += Tangent;
				VertexTangents[Index[2]] += Tangent;
			}

			// Normalize and set normals and tangents
			for (int32 i = 0; i < totalSize; i++)
			{
				VertexNormals[i].Normalize();
				//UE_LOG(LogTemp, Warning, TEXT("Normal %d: %s"), i, *VertexNormals[i].ToString());
				Builder.SetNormal(i, VertexNormals[i]);

				VertexTangents[i].Normalize();
				Builder.SetTangent(i, VertexTangents[i]);
			}
			
			Component->SetIsActive(true);
			if (!Component->GetIsInitialised()) {
				FString ComponentName = FString::Printf(TEXT("MainSection %i"), Component->GetIndex());
				FString ComponentName2 = FString::Printf(TEXT("MeshSection %i"), Component->GetIndex());
				FRealtimeMeshLODKey keyLOD = FRealtimeMeshLODKey::FRealtimeMeshLODKey(0);
				FRealtimeMeshSectionGroupKey GroupKey = FRealtimeMeshSectionGroupKey::Create(keyLOD, FName(ComponentName));
				RealtimeMesh->CreateSectionGroup(GroupKey, *StreamSet);
				const FRealtimeMeshSectionKey Key = FRealtimeMeshSectionKey::Create(GroupKey, FName(ComponentName2));
				Component->SetGroupKey(GroupKey);
				Component->SetKey(Key);
				const FRealtimeMeshStreamRange StreamRange(0, totalSize - 1, 0, totalSize - 1);
				FRealtimeMeshSectionConfig sectionCongig;
				sectionCongig.bCastsShadow = false;
				sectionCongig.MaterialSlot = 0;
				sectionCongig.DrawType = ERealtimeMeshSectionDrawType::Static;
				bool hasCollision = Component->GetLOD() == 0;
				RealtimeMesh->CreateSection(Key, sectionCongig, StreamRange, hasCollision);
				RealtimeMesh->UpdateSectionConfig(Key, sectionCongig, hasCollision);
				RealtimeMesh->UpdateSectionGroup(GroupKey, *StreamSet);
				Component->SetIsInitialised(true);
			}
			else if (RealtimeMesh != nullptr) {
				FRealtimeMeshSectionConfig sectionCongig;
				sectionCongig.bCastsShadow = false;
				sectionCongig.MaterialSlot = 0;
				sectionCongig.DrawType = ERealtimeMeshSectionDrawType::Static;
				bool hasCollision = Component->GetLOD() == 0;
				RealtimeMesh->UpdateSectionConfig(Component->GetKey(), sectionCongig, hasCollision);
				RealtimeMesh->UpdateSectionGroup(Component->GetGroupKey(), *StreamSet);
			}
			int32 VertexCount = Builder.NumVertices();
			//UE_LOG(LogTemp, Display, TEXT("SET VERTS: %i, SET LOD %i"), VertexCount, Component->GetLOD());
			Component->SetOldLOD(Component->GetLOD());
			Component->SetPrevBuilder(BuilderPtr);
			Component->SetPrevStreamset(MoveTemp(StreamSet));
		});
}

void AProceduralTerrain::GeneratePathMesh(int32 offset)
{

	////////////////////////////////////
	TSharedPtr<TRealtimeMeshBuilderLocal<uint16, FPackedNormal, FVector2DHalf, 1>> StoredBuilder = (offset==0)?nullptr:Path->GetPathBuilder();
	TSharedPtr<FRealtimeMeshStreamSet> StreamSet = (offset == 0) ? MakeShared<FRealtimeMeshStreamSet>(): Path->GetPathStreamSet();
	auto BuilderPtr = MakeShared<TRealtimeMeshBuilderLocal<uint16, FPackedNormal, FVector2DHalf, 1>>(*StreamSet);
	auto& Builder = (offset==0) ? *BuilderPtr : *StoredBuilder.Get();

	Builder.EnableTangents();
	Builder.EnableTexCoords();
	Builder.EnablePolyGroups();

	for (int32 i = offset; i < (PathPoints.Num() - 2); ++i)
	{

		FVector3f StartPoint = PathPoints[i];
		FVector3f EndPoint = PathPoints[(i + 1)];
		float StartHeightCenter = StartPoint.Z;
		float EndHeightCenter = EndPoint.Z;
		FVector3f CurrentPoint3D = FVector3f(StartPoint.X, StartPoint.Y, StartHeightCenter + HeightAdjust);
		FVector3f NextPoint3D = FVector3f(EndPoint.X, EndPoint.Y, EndHeightCenter + HeightAdjust);
		FVector3f Forward = (NextPoint3D - CurrentPoint3D).GetSafeNormal();
		FVector3f Up = FVector3f::UpVector;
		FVector3f Right = FVector3f::CrossProduct(Up, Forward).GetSafeNormal();

		TArray<FVector3f> CurrentNextVertexes;
		int32 IndexOffset = Builder.NumVertices();
		for (int l = 0; l < (ThicknessDetail * 2) + 1; l++) {
			float UCoord = static_cast<float>(l) / (ThicknessDetail * 2);
			float VCoord = i * PathTextureScale;
			if (l == ThicknessDetail) {
				if (i == 0) {
					Builder.AddVertex(NextPoint3D).SetTexCoord(FVector2f(UCoord, VCoord));
					Builder.AddVertex(CurrentPoint3D).SetTexCoord(FVector2f(UCoord, (i + 1) * PathTextureScale));
					CurrentNextVertexes.Add(NextPoint3D);
				}
				else {
					Builder.AddVertex(NextPoint3D).SetTexCoord(FVector2f(UCoord, VCoord));
					Builder.AddVertex(LastVertexes[l]).SetTexCoord(FVector2f(UCoord, (i + 1) * PathTextureScale));
					CurrentNextVertexes.Add(NextPoint3D);
				}
			}
			else if (l < ThicknessDetail) {
				// left
				float EndHeightLeft = CalculateNoiseAtPoint((EndPoint.X - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.X) / Scale, (EndPoint.Y - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.Y) / Scale);
				EndHeightLeft = FMath::Lerp(EndHeightLeft, EndHeightCenter, Flatness);
				FVector3f EndLeft = FVector3f(EndPoint.X - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.X, EndPoint.Y - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.Y, EndHeightLeft + HeightAdjust);
				CurrentNextVertexes.Add(EndLeft);
				if (i == 0) {
					float StartHeightLeft = CalculateNoiseAtPoint((StartPoint.X - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.X) / Scale, (StartPoint.Y - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.Y) / Scale);
					StartHeightLeft = FMath::Lerp(StartHeightLeft, StartHeightCenter, Flatness);
					FVector3f StartLeft = FVector3f(StartPoint.X - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.X, StartPoint.Y - ((Thickness / ThicknessDetail) * (ThicknessDetail - l)) * Right.Y, StartHeightLeft + HeightAdjust);

					Builder.AddVertex(EndLeft).SetTexCoord(FVector2f(UCoord, VCoord));
					Builder.AddVertex(StartLeft).SetTexCoord(FVector2f(UCoord, (i + 1) * PathTextureScale));
				}
				else {
					Builder.AddVertex(EndLeft).SetTexCoord(FVector2f(UCoord, VCoord));
					Builder.AddVertex(LastVertexes[l]).SetTexCoord(FVector2f(UCoord, (i + 1) * PathTextureScale));
				}
			}
			else {
				//right
				float EndHeightRight = CalculateNoiseAtPoint((EndPoint.X + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.X) / Scale, (EndPoint.Y + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.Y) / Scale);
				EndHeightRight = FMath::Lerp(EndHeightRight, EndHeightCenter, Flatness);
				FVector3f EndRight = FVector3f(EndPoint.X + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.X, EndPoint.Y + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.Y, EndHeightRight + HeightAdjust);
				CurrentNextVertexes.Add(EndRight);
				if (i == 0) {

					float StartHeightRight = CalculateNoiseAtPoint((StartPoint.X + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.X) / Scale, (StartPoint.Y + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.Y) / Scale);
					StartHeightRight = FMath::Lerp(StartHeightRight, StartHeightCenter, Flatness);
					FVector3f StartRight = FVector3f(StartPoint.X + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.X, StartPoint.Y + ((Thickness / ThicknessDetail) * (l - ThicknessDetail)) * Right.Y, StartHeightRight + HeightAdjust);

					Builder.AddVertex(EndRight).SetTexCoord(FVector2f(UCoord, VCoord));
					Builder.AddVertex(StartRight).SetTexCoord(FVector2f(UCoord, (i + 1) * PathTextureScale));
				}
				else {
					Builder.AddVertex(EndRight).SetTexCoord(FVector2f(UCoord, VCoord));
					Builder.AddVertex(LastVertexes[l]).SetTexCoord(FVector2f(UCoord, (i + 1) * PathTextureScale));
				}
			}
		}
		LastVertexes = CurrentNextVertexes;

		for (int j = 0; j < (ThicknessDetail * 2); j++) {
			int32 point1 = (IndexOffset + 1 + (j * 2));
			int32 point2 = (IndexOffset + 2 + (j * 2));
			int32 point3 =(IndexOffset + 0 + (j * 2));

			int32 point4 = (IndexOffset + 3 + (j * 2));
			int32 point5 = (IndexOffset + 2 + (j * 2));
			int32 point6 = (IndexOffset + 1 + (j * 2));

			Builder.AddTriangle(point1, point2, point3);
			Builder.AddTriangle(point4, point5, point6);
		}

	}


	for (int i = prevEnd; i < Builder.NumVertices(); i) {
		FVector3f getPos = Builder.GetPosition(i);
		getPos.Z += EdgeHeightOffset;
		Builder.SetPosition(i, getPos);
		FVector3f getPos2 = Builder.GetPosition(i+1);
		getPos2.Z += EdgeHeightOffset;
		Builder.SetPosition(i+1, getPos2);
		FVector3f getPosThick = Builder.GetPosition(i + (((ThicknessDetail * 2) + 1) * 2) - 2);
		getPosThick.Z += EdgeHeightOffset;
		Builder.SetPosition(i + (((ThicknessDetail * 2) + 1) * 2) - 2, getPosThick);
		FVector3f getPosThick2 = Builder.GetPosition(i + (((ThicknessDetail * 2) + 1) * 2) - 1);
		getPosThick2.Z += EdgeHeightOffset;
		Builder.SetPosition(i + (((ThicknessDetail * 2) + 1) * 2) - 1, getPosThick2);

		i += ((ThicknessDetail * 2) + 1) * 2;
	}

	int VertsPerPoint = (((ThicknessDetail * 2) + 1) * 2);
	

	for (int i = prevEnd; i < Builder.NumVertices(); i++) {
		int32 GridX = FMath::FloorToInt(PathPoints[pointCount].X / currentGridSize);
		int32 GridY = FMath::FloorToInt(PathPoints[pointCount].Y / currentGridSize);
		//UE_LOG(LogTemp, Display, TEXT("Point at Grid (%i,%i)"), GridX, GridY);
		FIntPoint Cell = FIntPoint(GridX, GridY);
		PathVertMap.FindOrAdd(Cell).Add(Builder.GetPosition(i));
		vertsInPoint++;
		if (vertsInPoint >= VertsPerPoint) {
			pointCount++;
			vertsInPoint = 0;
		}
	}
	pointCount+=2;
	UE_LOG(LogTemp, Display, TEXT("Points %i, Verts %i, Count %i"), PathPoints.Num(), Builder.NumVertices(), pointCount);
	//pointCount++;
	prevEnd = Builder.NumVertices();


	if (offset == 0) {
		Path = new PathComponent(nullptr);
		FString ComponentName = FString::Printf(TEXT("Path"));
		FString ComponentName2 = FString::Printf(TEXT("PathSection"));
		FRealtimeMeshLODKey keyLOD = FRealtimeMeshLODKey::FRealtimeMeshLODKey(0);
		FRealtimeMeshSectionGroupKey GroupKey = FRealtimeMeshSectionGroupKey::Create(keyLOD, FName(ComponentName));
		RealtimeMesh->CreateSectionGroup(GroupKey, *StreamSet);
		const FRealtimeMeshSectionKey Key = FRealtimeMeshSectionKey::Create(GroupKey, FName(ComponentName2));
		Path->SetGroupKey(GroupKey);
		Path->SetKey(Key);
		const FRealtimeMeshStreamRange StreamRange(0, Builder.NumVertices() - 1, 0, Builder.NumVertices() - 1);
		FRealtimeMeshSectionConfig sectionCongig;
		sectionCongig.bCastsShadow = false;
		sectionCongig.MaterialSlot = 1;
		sectionCongig.DrawType = ERealtimeMeshSectionDrawType::Static;
		RealtimeMesh->CreateSection(Key, sectionCongig, StreamRange, true);
		RealtimeMesh->UpdateSectionGroup(GroupKey, *StreamSet);
		Path->SetPathBuilder(BuilderPtr);
		Path->SetPathStreamset(MoveTemp(StreamSet));
	}
	else {
		FRealtimeMeshSectionConfig sectionCongig;
		sectionCongig.bCastsShadow = false;
		sectionCongig.MaterialSlot = 1;
		sectionCongig.DrawType = ERealtimeMeshSectionDrawType::Static;
		sectionCongig.bForceOpaque = true;
		RealtimeMesh->UpdateSectionConfig(Path->GetKey(), sectionCongig, true);
		RealtimeMesh->UpdateSectionGroup(Path->GetGroupKey(), *StreamSet);
	}


}
